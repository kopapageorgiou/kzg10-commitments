from typing import List, NamedTuple, Tuple, Union
from random import randint
from BlockchainLib import SmartContract
from BLS import *
from py_ecc.fields import bn128_FQ2 as FQ2
from py_ecc.typing import Point2D, Field
from KZG10 import *
import time

def test_single_proof_lagrange():
	values = [randint(1,50) for _ in range(10)]
	index = 2
	kzg = KZG10()
	smartContract = SmartContract()
	lagrange = kzg.choose_method(KZG10.LAGRANGE)
	#start = time.process_time()
	coeffs = lagrange.interpolate(values)
	commitment = lagrange.generate_commitment(coeffs)
	expected_commitment = smartContract.commit([int(c) for c in coeffs])
	assert commitment == expected_commitment, "Commitment is not correct"
	#end = time.process_time() - start
	proof = lagrange.generate_proof(coeffs, index)
	value_y = lagrange.eval_poly_at(coeffs, index)
	expected_value_y = smartContract.evalPolyAt([int(c) for c in coeffs], index)
	assert value_y == expected_value_y, "Evaluation is not correct"
	verification = kzg.verify_off_chain(commitment, proof, index, value_y)
	expected_verification = smartContract.verify(format_FQ_G1Point(commitment), format_FQ_G1Point(proof), index, format_field_to_int(value_y))
	assert verification == expected_verification, "Verification is not correct"
	assert verification == True, "Verification is not correct"
	#print(f"Lagrange Interpolation for {len(values)} values")
	#print('-'*50)
	#print(f"Cpu time for commitment generation: {end} secs")
	#print(f"Is evaluation correct? {int(value_y) == values[index]}")
	#print(f"Is pairing correct? {verif}")
	#print('='*50)

def testMonomial(x_values, values, index):
	kzg = KZG10()
	monomial = kzg.choose_method(KZG10.MONOMIAL)
	start = time.process_time()
	coeffs = monomial.interpolate(values)
	print (f"coeffs {coeffs}")
	
	#print("y:", value_y)
	commit = monomial.generate_commitment(coeffs)
	end = time.process_time() - start
	proof = monomial.generate_proof(coeffs, index)
	value_y = monomial.eval_poly_at(coeffs, index)	
	verif = kzg.verify_off_chain(commit, proof, index, value_y)
	print(f"Monomial Interpolation for {len(values)} values")
	print('-'*50)
	print(f"Cpu time for commitment generation: {end} secs")
	print(f"Is evaluation correct? {int(value_y) == values[index]}")
	print(f"Is pairing correct? {verif}")
	print('='*50)

def test_single_proof_newtons():
	values = [randint(1,50) for _ in range(10)]
	index = 2
	kzg = KZG10()
	smartContract = SmartContract()
	newton = kzg.choose_method(KZG10.NEWTON)
	#print("values=",values)
	#start = time.process_time()
	coeffs = newton.interpolate(values)
	#print("coeffs:", coeffs)
	
	#print("y:", y)
	commitment = newton.generate_commitment(coeffs)
	expected_commitment = smartContract.commit([int(c) for c in coeffs])
	assert commitment == expected_commitment, "Commitment is not correct"
	#end = time.process_time() - start
	proof = newton.generate_proof(coeffs, index)
	value_y = newton.eval_poly_at(coeffs, index)
	expected_value_y = smartContract.evalPolyAt([int(c) for c in coeffs], index)
	assert value_y == expected_value_y, "Evaluation is not correct"
	#print (f"commit : {commit},  proof : {proof}, index: {index} , value_y : {y}")
	
	verification = kzg.verify_off_chain(commitment, proof, index, value_y)
	expected_verification = smartContract.verify(format_FQ_G1Point(commitment), format_FQ_G1Point(proof), index, format_field_to_int(value_y))
	assert verification == expected_verification, "Verification is not correct"
	# print(f"Newton Interpolation for {len(values)} values")
	# print('-'*50)
	# print(f"Cpu time for commitment generation: {end} secs")
	# print(f"Is evaluation correct? {int(y) == values[index]}")
	# print(f"Is pairing correct? {verification}")
	# print('='*50)

def test_multi_proof_newtons_with_correct_indices():
	values = [randint(1,50) for _ in range(10)]
	indices = [0,1,2,3]
	kzg = KZG10()
	smartContract = SmartContract()
	newton = kzg.choose_method(KZG10.NEWTON)
	coeffs = newton.interpolate(values)
	#print("coeffs: ", coeffs)
	commitment = newton.generate_commitment(coeffs)
	expected_commitment = smartContract.commit([int(c) for c in coeffs])
	assert commitment == expected_commitment, "Commitment is not correct"
	for index in indices:
		value_y = newton.eval_poly_at(coeffs, index)
		expected_value_y = smartContract.evalPolyAt([int(c) for c in coeffs], index)
		assert value_y == expected_value_y, "Evaluation is not correct"
	
	multiproof, icoeff, zpoly = newton.generate_multi_proof(coeffs, indices, [values[i] for i in indices])
	#multiproof, icoeff, zpoly = kzg.generate_multi_proof(coeffs, indices, [values[i] for i in indices])
	#print("proof: ", multiproof)
	result = smartContract.verifyMulti(format_FQ_G1Point(commitment), format_proof(multiproof), indices, [values[i] for i in indices], [int(c) for c in icoeff], [int(z) for z in zpoly])
	assert result==True, "Verification is not correct"

def test_multi_proof_newtons_with_incorrect_indices_in_proof():
	values = [randint(1,50) for _ in range(10)]
	indices = [0,1,2,3]
	kzg = KZG10()
	smartContract = SmartContract()
	newton = kzg.choose_method(KZG10.NEWTON)
	coeffs = newton.interpolate(values)
	#print("coeffs: ", coeffs)
	commitment = newton.generate_commitment(coeffs)
	expected_commitment = smartContract.commit([int(c) for c in coeffs])
	assert commitment == expected_commitment, "Commitment is not correct"
	for index in indices:
		value_y = newton.eval_poly_at(coeffs, index)
		expected_value_y = smartContract.evalPolyAt([int(c) for c in coeffs], index)
		assert value_y == expected_value_y, "Evaluation is not correct"
	
	multiproof, icoeff, zpoly = newton.generate_multi_proof(coeffs, [2, 3, 4, 5], [values[i] for i in indices])
	try:
		result = smartContract.verifyMulti(format_FQ_G1Point(commitment), format_proof(multiproof), indices, [values[i] for i in indices], [int(c) for c in icoeff], [int(z) for z in zpoly])
		assert result==False, "Verification is not correct"
	except ValueError as e:
		assert e.args[0]['message'].__contains__("Verifier.verifyMulti: invalid _zCoeffs"), "Verification is not the expected result"

def test_multi_proof_newtons_with_incorrect_values_in_proof():
	values = [randint(1,50) for _ in range(10)]
	indices = [0,1,2,3]
	kzg = KZG10()
	smartContract = SmartContract()
	newton = kzg.choose_method(KZG10.NEWTON)
	coeffs = newton.interpolate(values)
	#print("coeffs: ", coeffs)
	commitment = newton.generate_commitment(coeffs)
	expected_commitment = smartContract.commit([int(c) for c in coeffs])
	assert commitment == expected_commitment, "Commitment is not correct"
	for index in indices:
		value_y = newton.eval_poly_at(coeffs, index)
		expected_value_y = smartContract.evalPolyAt([int(c) for c in coeffs], index)
		assert value_y == expected_value_y, "Evaluation is not correct"
	
	multiproof, icoeff, zpoly = newton.generate_multi_proof(coeffs, indices, [values[i] for i in [2, 3, 4, 5]])
	#multiproof, icoeff, zpoly = kzg.generate_multi_proof(coeffs, indices, [values[i] for i in indices])
	#print("proof: ", multiproof)
	result = smartContract.verifyMulti(format_FQ_G1Point(commitment), format_proof(multiproof), indices, [values[i] for i in indices], [int(c) for c in icoeff], [int(z) for z in zpoly])
	try:
		result = smartContract.verifyMulti(format_FQ_G1Point(commitment), format_proof(multiproof), indices, [values[i] for i in indices], [int(c) for c in icoeff], [int(z) for z in zpoly])
		assert result==False, "Verification is not correct"
	except ValueError:
		pass
	

def testSpline(x_values, values, index):
	kzg = KZG10()
	
	start = time.process_time()
	coeffs = kzg.b_spline_interpolation(x_values, values,0)
	#coeffs2 = kzg.pre_coeffs(coeffs, index)
	print("Values", values)
	print("Coeffs", coeffs, len(coeffs))
	y = kzg.evaluate_spline_point(coeffs, index)
	#ofcomt = kzg.generate_commitment(coeffs)
	comt = kzg.generate_commitment(coeffs)
	prf = kzg.custom_generate_proof(coeffs, index)
	end = time.process_time() - start
	print(f"Cubic Spline Interpolation for {len(values)} values")
	print('-'*50)
	#print(ofcomt)
	print(comt)
	print(prf)
	print("Y", y)
	print(f"Cpu time: {end} secs")
	print(f"Is evaluation correct? {int(y) == values[index]}")
	print(f"Is pairing correct? {kzg.verify_off_chain(comt, prf, index, y)}")

def testLinear(x_values, values, index):
	kzg = KZG10()
	start = time.process_time()
	coeffs = kzg.linear_interpolation(x_values, values)
	print("size" + str(len(coeffs)))
	y = kzg.linear_inter_evaluation(coeffs, x_values, index)
	#print(coefficients)
	end = time.process_time() - start
	print(f"Linear Interpolation for {len(values)} values")
	print('-'*50)
	print(f"Cpu time: {end} secs")
	print(f"Is evaluation correct? {int(y) == values[index]}")
	print(f"Is pairing correct? {kzg.verify_off_chain(comt, prf, index, values[index])}")

def testMonomial(x_values, values, index):
	kzg = KZG10()
	mon = kzg.choose_method(KZG10.MONOMIAL)
	start = time.process_time()
	coeffs = mon.interpolate(values)
	#print("coeffs:", coeffs)
	commit = mon.generate_commitment(coeffs)
	end = time.process_time() - start
	proof = mon.generate_proof(coeffs, index)
	value_y = mon.eval_monomial_at(coeffs, index)	
	print (f"commit : {commit},  proof : {proof}, index: {index} , value_y : {value_y}")
	verif = kzg.verify_off_chain(commit, proof, index, value_y)
	print(f"Monomial Interpolation for {len(values)} values")
	print('-'*50)
	print(f"Cpu time for commitment generation: {end} secs")
	print(f"Is evaluation correct? {int(value_y) == values[index]}")
	print(f"Is pairing correct? {verif}")
	print('='*50)


if __name__ == "__main__":
	#testRun()
	#testMultiProof()
	#values = [randint(1,50) for i in range(10)]
	#values = [1,5,17]
	
	#x_values = [i for i in range(len(values))]
	#x_values = [0,1,2]
	#testLagrange(x_values, values, 2)
	#testNewton(values,2)
	#testMonomial(x_values,values, 2)
	#testMultiProofNewton(x_values,values, [0,1,2,3])
	test_multi_proof_newtons_with_incorrect_indices_in_proof()
	#testSpline(x_values, values, 0)
	#testLinear(x_values, values, randint(0, len(values)-1))
