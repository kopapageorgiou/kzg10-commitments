from py_ecc.bn128.bn128_curve import multiply, G1, G2, add
from py_ecc.bls import G2ProofOfPossession as bls_pop
from py_ecc.bn128 import pairing
import os
from web3 import Web3
#from smartContractInteraction import smartContractInteraction

"""
#! Method to create BLS private and public key pair
#* @param fixed: fixed integer for private key
"""
def GenerateKeys(fixed = None):
    privateKey = fixed or bls_pop.KeyGen(os.urandom(32))
    publicKey = getPublicKey(privateKey)

    return (privateKey, publicKey)

"""
#! Method to get the point of public key on bn128 curve from private key
#* @param private_key: The BLS private key
"""
def getPublicKey(private_key: int):
    return multiply(G1, private_key)

"""
#! Method to sign a message with BLS private key 
#* @param names: The signers
#* @param status: The message to be signed
#* @param private_key: the private key to sign the message
#* @return: The message hashed on curve (G2 Point), the message hashed (bytes)
"""
def signBLS(address1: str, funcName: str, arg: int, privateKey: int):
    #addr1 = Web3.toChecksumAddress(address1)
    addr = address1.lower()
    conc = Web3.solidityKeccak(['string','string', 'string'], [addr, funcName, str(arg)])

    data = hashMessageOnCurve(_bytesToInt(conc))
    return multiply(data, privateKey), conc

"""
#! Method to return the point of the message on the bn128 curve
#* @param public keys: The BLS public keys list
"""
def hashMessageOnCurve(message: int):
    return multiply(G2, message)

"""
#! Method to aggregate signatures 
#* @param signatures: The BLS signatures list
"""
def getAggregatedSignature(sigs: list):
    aggregated = add(sigs[0], sigs[1])
    for i in range(2, len(sigs)):
        aggregated = add(aggregated, sigs[i])
    print(aggregated)
    return formatG2(aggregated)


"""
#! Method to aggregate public keys 
#* @param public keys: The BLS public keys list
"""
def getAggregatedPublicKey(public_keys: list):
    assert(len(public_keys)>=2)
    aggregatedPublicKey = add(public_keys[0], public_keys[1])
    for i in range(2,len(public_keys)):
        aggregatedPublicKey = add(aggregatedPublicKey, public_keys[i])
    
    return formatG2(aggregatedPublicKey)

"""
#! Method to convert bytes to decimal 
#* @param data: The message in bytes format
"""
def _bytesToInt(data: bytes):
    return int.from_bytes(data, "big")
    
"""
#! Method to verify aggregated signatures calling method from smart contract
#* @param smart_contract: The instance of smart contract
#* @param contract address to be called on-chain
#* @param function name to be called
#* @param arg the argument of the funcion
#* @param address2: The address of the caller
#* @param status: The hash message to be signed
#* @param aggregated_signature: The BLS aggregated signature
#* @param public_keys: The list of public keys
"""
#def verifyOnChain(smart_contract: smartContractInteraction, address: str, funcName: str, arg: int, address2: str, status: bytes, aggregated_signature: tuple, public_key: list[tuple]):
#    addr = Web3.toChecksumAddress(address)
#    addr2 = Web3.toChecksumAddress(address2)
#    return smart_contract.verify(addr, funcName, arg, addr2, status, aggregated_signature, public_key)

"""
#! Method to verify aggregated signatures N of M public keys calling method from smart contract
#? This method is under development
#* @param smart_contract: The instance of smart contract
#* @param bitmask: The number to be used as mask for public key selection
#* @param signers: The signers list
#* @param status: The message signed
#* @param signatures: The list of signatures
"""
#def verifyNofMOnChain(smart_contract: smartContractInteraction, bitmask: int, signers: list, status: bytes, signatures: list):
#    data = _bytesToInt(status)
#    return smart_contract.verify2NofM(bitmask, signers, data, signatures)

"""
#! Method to verify aggregated signatures without using smart contract
#* @param aggregated_public_key: The BLS aggregated public key
#* @param status: The message signed
#* @param aggregated_signature: The BLS aggregated signature
"""
def verifyOffChain(aggregated_public_key: tuple, status: bytes, aggregated_signature: tuple):
    data = hashMessageOnCurve(_bytesToInt(status))
    pairing1 = pairing(data, aggregated_public_key)
    pairing2 = pairing(aggregated_signature, G1)
    return pairing1 == pairing2

"""
#! Method to format a G2 point in order to be recognized from the smart contract's methods
#* @param data: The G2 point
"""
def formatG2(data):
    p1, p2 = data
    x1, x2 = p1.coeffs
    y1, y2 = p2.coeffs
    return ([int(x1),int(x2)],[int(y1),int(y2)])

"""
#! Method to format a G1 point in order to be recognized from the smart contract's methods
#* @param data: The G1 point
"""
def formatG1(data):
    x, y = data
    return (int(x), int(y))

def formatG1_FQ(data):
    x, y = data
    print((int(str(x)), int(str(y))))
    return (int(str(x)), int(str(y)))