from web3 import Web3
import json, configparser
class SmartContract(object):
    
    def __init__(self, settings='settings.ini') -> None:
        
        config = configparser.ConfigParser()
        try:
            config.read(settings)
            addr = f"{config['SERVER']['Address']}:{config['SERVER']['Port']}"
            timeout = int(config['SERVER']['Timeout'])
            path = config['SMART CONTRACT']["Abi"]
        except:
            print('An error occured while reading config file')
        try:
            self.web3 = Web3(Web3.HTTPProvider(addr, request_kwargs={'timeout':timeout}))
            address = self.web3.to_checksum_address(config['SMART CONTRACT']['Address'])
        except:
            print("Could not establish connection with the Blockchain test network")
        try:
            with open("Blockchain/.build/Verifier.json", 'r') as fp:
                json_info = json.load(fp)
                abi = json_info['abi']
                self.account = self.web3.eth.accounts[0]
                self.web3.eth.default_account = self.account

                self.contract = self.web3.eth.contract(address=address, abi=abi)
        except Exception as e:
            print('An error occured while loading abi from file',e )

    def verify(self, commitment, proof, index, value):
        return self.contract.functions.verify(commitment, proof, index, value).call()
    
    def verifyMulti(self, commitment, proof, indices, values, iCoeffs, zCoeffs):
        return self.contract.functions.verifyMulti(commitment, proof, indices, values, iCoeffs, zCoeffs).call()

    def commit(self, coefficients):
        return self.contract.functions.commit(coefficients).call()

    def evalPolyAt(self, coefficients, index):
        return self.contract.functions.evalPolyAt(coefficients, index).call()
    
    def evalPolyAt2(self, coefficients, index):
        return self.contract.functions.evalPolyAt2(coefficients, index).call()
    
    def pairing(self, a1, a2, b1, b2):
        return self.contract.functions.pairingTest(a1,a2,b1,b2).call()
    
    def testAdd(self, coeffs):
        return self.contract.functions.testAdd(coeffs).call()

    def testMul(self, s):
        return self.contract.functions.testMul(s).call()