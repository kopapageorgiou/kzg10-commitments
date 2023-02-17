// Modified from https://github.com/appliedzkp/semaphore/blob/master/contracts/sol/verifier.sol
pragma experimental ABIEncoderV2;
pragma solidity ^0.8.0;

//import "./Pairing.sol";
import "./BLSVerification.sol";
import { Constants } from "./Constants.sol";

contract Verifier is Constants {
    uint256 private constant FIELD_MODULUS = 0x30644e72e131a029b85045b68181585d97816a916871ca8d3c208c16d87cfd47;
    //using Pairing for *;
    //using BLSVerification for *;

    // The G1 generator
    BLSVerification.G1Point SRS_G1_0 = BLSVerification.G1Point({
        X: Constants.SRS_G1_X[0],
        Y: Constants.SRS_G1_Y[0]
    });

    // The G2 generator
    BLSVerification.G2Point g2Generator = BLSVerification.G2Point({
        X: [ Constants.SRS_G2_X_0[0], Constants.SRS_G2_X_1[0] ],
        Y: [ Constants.SRS_G2_Y_0[0], Constants.SRS_G2_Y_1[0] ]

    });

    BLSVerification.G2Point SRS_G2_1 = BLSVerification.G2Point({
        X: [ Constants.SRS_G2_X_0[1], Constants.SRS_G2_X_1[1] ],
        Y: [ Constants.SRS_G2_Y_0[1], Constants.SRS_G2_Y_1[1] ]
    });

    function pairingTest(
        BLSVerification.G1Point memory a1,
        BLSVerification.G2Point memory a2,
        BLSVerification.G1Point memory b1,
        BLSVerification.G2Point memory b2
    ) public returns (bool) {
        return BLSVerification.pairing2(a1, a2, b1, b2);
    }

    /*
     * Verifies a single-point evaluation of a polynominal using the KZG
     * commitment scheme.
     *    - p(X) is a polynominal
     *    - _value = p(_index) 
     *    - commitment = commit(p)
     *    - proof = genProof(p, _index, _value)
     * Returns true if and only if the following holds, and returns false
     * otherwise:
     *     e(commitment - commit([_value]), G2.g) == e(proof, commit([0, 1]) - zCommit)
     * @param _commitment The KZG polynominal commitment.
     * @param _proof The proof.
     * @param _index The x-value at which to evaluate the polynominal.
     * @param _value The result of the polynominal evaluation.
     */
    function verify(
        BLSVerification.G1Point memory _commitment,
        BLSVerification.G1Point memory _proof,
        uint256 _index,
        uint256 _value
    ) public returns (bool) {
        // Make sure each parameter is less than the prime q
        require(_commitment.X < BABYJUB_P, "Verifier.verifyKZG: _commitment.X is out of range");
        require(_commitment.Y < BABYJUB_P, "Verifier.verifyKZG: _commitment.Y is out of range");
        require(_proof.X < BABYJUB_P, "Verifier.verifyKZG: _proof.X is out of range");
        require(_proof.Y < BABYJUB_P, "Verifier.verifyKZG: _proof.Y is out of range");
        require(_index < BABYJUB_P, "Verifier.verifyKZG: _index is out of range");
        require(_value < BABYJUB_P, "Verifier.verifyKZG: _value is out of range");

        // Check that 
        //     e(commitment - aCommit, G2.g) == e(proof, xCommit - zCommit)
        //     e(commitment - aCommit, G2.g) / e(proof, xCommit - zCommit) == 1
        //     e(commitment - aCommit, G2.g) * e(proof, xCommit - zCommit) ^ -1 == 1
        //     e(commitment - aCommit, G2.g) * e(-proof, xCommit - zCommit) == 1
        // where:
        //     aCommit = commit([_value]) = SRS_G1_0 * _value
        //     xCommit = commit([0, 1]) = SRS_G2_1
        //     zCommit = commit([_index]) = SRS_G2_1 * _index

        // To avoid having to perform an expensive operation in G2 to compute
        // xCommit - zCommit, we instead check the equivalent equation:
        //     e(commitment - aCommit, G2.g) * e(-proof, xCommit) * e(-proof, -zCommit) == 1
        //     e(commitment - aCommit, G2.g) * e(-proof, xCommit) * e(proof, zCommit) == 1
        //     e(commitment - aCommit, G2.g) * e(-proof, xCommit) * e(index * proof, G2.g) == 1
        //     e((index * proof) + (commitment - aCommit), G2.g) * e(-proof, xCommit) == 1


        // Compute commitment - aCommitment
        BLSVerification.G1Point memory commitmentMinusA = BLSVerification.add(
            _commitment,
            BLSVerification.negate(
                BLSVerification.mul(SRS_G1_0, _value)
            )
        );

        // Negate the proof
        BLSVerification.G1Point memory negProof = BLSVerification.negate(_proof);

        // Compute index * proof
        BLSVerification.G1Point memory indexMulProof = BLSVerification.mul(_proof, _index);

        // Returns true if and only if
        // e((index * proof) + (commitment - aCommitment), G2.g) * e(-proof, xCommit) == 1
        return BLSVerification.pairing2(
            BLSVerification.add(indexMulProof, commitmentMinusA),
            g2Generator,
            negProof,
            SRS_G2_1
        );
    }

    function evalPolyAt2(uint256[] memory coeffs, uint256 x) public view returns (uint256) {
        uint256 result = 0;
        uint256 n = coeffs.length;
        uint256 fieldSize = Constants.BABYJUB_P;
        for (uint256 i = 0; i < n; i++) {
            uint256 term = mulmod(coeffs[i], pow(x, n-i-1, fieldSize), fieldSize);
            result = addmod(result, term, fieldSize);
        }
        
        return result;
    }

    function pow(uint256 base, uint256 exp, uint256 fieldSize) internal view returns (uint256) {
        if (exp == 0) {
            return 1;
        } else if (exp == 1) {
            return base;
        } else {
            uint256 temp = pow(base, exp / 2, fieldSize);
            if (exp % 2 == 0) {
                return mulmod(temp, temp, fieldSize);
            } else {
                return mulmod(mulmod(temp, temp, fieldSize), base, fieldSize);
            }
        }
    }

    /*
     * @return A KZG commitment to a polynominal
     * @param coefficients The coefficients of the polynomial to which to
     *                     commit.
     */
    function commit(
        uint256[] memory coefficients
    ) public returns (BLSVerification.G1Point memory) {

        BLSVerification.G1Point memory result = BLSVerification.G1Point(0, 0);

        for (uint256 i = 0; i < coefficients.length; i ++) {
            result = BLSVerification.add(
                result,
                BLSVerification.mul(
                    BLSVerification.G1Point({
                        X: Constants.SRS_G1_X[i],
                        Y: Constants.SRS_G1_Y[i]
                    }),
                    coefficients[i]
                )
            );
        }
        return result;
    }
    function testAdd(
        uint256[] memory coefficients
    ) public returns (BLSVerification.G1Point memory) {

        BLSVerification.G1Point memory result = BLSVerification.G1Point(0, 0);

        for (uint256 i = 0; i < coefficients.length; i ++) {
            result = BLSVerification.add(
                result,
                    BLSVerification.G1Point({
                        X: Constants.SRS_G1_X[i],
                        Y: Constants.SRS_G1_Y[i]
                    })
            );
        }
        return result;
    }
    function testMul(uint256[] memory coeffs) public returns(BLSVerification.G1Point memory){
        BLSVerification.G1Point memory result = BLSVerification.G1Point(0, 0);
        for (uint256 i = 0; i < coeffs.length; i ++) {
            BLSVerification.G1Point memory temp = BLSVerification.G1Point(SRS_G1_X[i], SRS_G1_Y[i]);
            BLSVerification.G1Point memory mulTemp = BLSVerification.mul(temp, coeffs[i]);
            result = BLSVerification.add(result, mulTemp);
        }
        return result;
    }
    function _FQ2Add(
        BLSVerification.G1Point memory p1,
        BLSVerification.G1Point memory p2
    ) internal pure returns(BLSVerification.G1Point memory) {
        return (BLSVerification.G1Point({
            X: addmod(p1.X, p2.X, FIELD_MODULUS),
            Y: addmod(p1.Y, p2.Y, FIELD_MODULUS)
        })
        );
    }

    /*
     * @return The polynominal evaluation of a polynominal with the specified
     *         coefficients at the given index.
     */
    function evalPolyAt(
        uint256[] memory _coefficients,
        uint256 _index
    ) public pure returns (uint256) {

        uint256 m = Constants.BABYJUB_P;
        uint256 result = 0;
        uint256 powerOfX = 1;

        for (uint256 i = 0; i < _coefficients.length; i ++) {
            uint256 coeff = _coefficients[i];
            assembly {
                result:= addmod(result, mulmod(powerOfX, coeff, m), m)
                powerOfX := mulmod(powerOfX, _index, m)
            }
        }
        return result;
    }
    
    /*
     * Verifies the evaluation of multiple points of a polynominal using the
     * KZG commitment scheme.
     *    - p(X) is a polynominal
     *    - commitment = commit(p)
     *    - For each y in _values and each x in _indices, y = p(x)
     *    - proof = genMultiProof(p, _indices)
     * Returns true if and only if the following holds, and returns false
     * otherwise:
     *     e(-commit(zPoly), proof) * e(commitment - commit(iPoly), g) == 1
     * @param _commitment The polynominal commitment.
     * @param _proof The proof.
     * @param _indices The x-values at which to evaluate the polynominal.
     * @param _values The evaluation of the polynominal at each index.
     * @param _iCoeffs The coefficients of a polynominal which interpolates
     *                 each index and corresponding y-value.
     * @param _zCoeffs The coefficients of a polynominal which intersects y=0
     *                 for each index.
     */
    function verifyMulti(
        BLSVerification.G1Point memory _commitment,
        BLSVerification.G2Point memory _proof,
        uint256[] memory _indices,
        uint256[] memory _values,
        uint256[] memory _iCoeffs,
        uint256[] memory _zCoeffs
    ) public returns (bool) {
        // Perform range checks
        require(_commitment.X < BABYJUB_P, "Verifier.verifyMultiKZG: _commitment.X is out of range");
        require(_commitment.Y < BABYJUB_P, "Verifier.verifyMultiKZG: _commitment.Y is out of range");
        require(_proof.X[0] < BABYJUB_P, "Verifier.verifyKZG: _proof.X[0] is out of range");
        require(_proof.X[1] < BABYJUB_P, "Verifier.verifyKZG: _proof.X[1] is out of range");
        require(_proof.Y[0] < BABYJUB_P, "Verifier.verifyKZG: _proof.Y[0] is out of range");
        require(_proof.Y[1] < BABYJUB_P, "Verifier.verifyKZG: _proof.Y[1] is out of range");

        for (uint256 i = 0; i < _iCoeffs.length; i ++) {
            require(_iCoeffs[i] < BABYJUB_P, "Verifier.verifyKZG: an _iCoeffs value is out of range");
        }

        for (uint256 i = 0; i < _zCoeffs.length; i ++) {
            require(_zCoeffs[i] < BABYJUB_P, "Verifier.verifyKZG: an _zCoeffs value is out of range");
        }

        // Check whether _iCoeffs and _zCoeffs are valid
        for (uint256 i = 0; i < _indices.length; i ++) {
            uint256 index = _indices[i];
            uint256 value = _values[i];
            require(index < BABYJUB_P, "Verifier.verifyKZG: an index is out of range");
            require(value < BABYJUB_P, "Verifier.verifyKZG: a value is out of range");

            uint256 zEval = evalPolyAt(_zCoeffs, _indices[i]);
            require(zEval == 0, "Verifier.verifyMulti: invalid _zCoeffs");

            uint256 iEval = evalPolyAt(_iCoeffs, _indices[i]);
            require(iEval == _values[i], "Verifier.verifyMulti: invalid _iCoeffs");
        }

        // Generate the KZG commitments to the i and z polynominals
        BLSVerification.G1Point memory zCommit = commit(_zCoeffs);
        BLSVerification.G1Point memory iCommit = commit(_iCoeffs);

        // Compute commitment - commit(iPoly)
        BLSVerification.G1Point memory commitmentMinusICommit =
            BLSVerification.add(
                _commitment,
                BLSVerification.negate(iCommit)
            );

        // Perform the pairing check
        return BLSVerification.pairing2(
            BLSVerification.negate(zCommit),
            _proof,
            commitmentMinusICommit,
            g2Generator
        );
    }

    /*
    // Uncomment to perform gas benchmarks
    function commitBenchmark(
        uint256[] memory _coefficients
    ) public {
        commit(_coefficients);
    }

    function verifyMultiBenchmark(
        Pairing.G1Point memory _commitment,
        Pairing.G2Point memory _proof,
        uint256[] memory _indices,
        uint256[] memory _values,
        uint256[] memory _iCoeffs,
        uint256[] memory _zCoeffs
    ) public {
        verifyMulti(_commitment, _proof, _indices, _values, _iCoeffs, _zCoeffs);
    }
    */

    
}
