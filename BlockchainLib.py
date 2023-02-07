from web3 import Web3
from brownie import accounts
import json, configparser
class smartContract(object):
    
    def __init__(self, settings='settings.ini') -> None:
        
        config = configparser.ConfigParser()
        try:
            config.read(settings)
            addr = f"{config['SERVER']['Address']}:{config['SERVER']['Port']}"
            print(addr)
            timeout = int(config['SERVER']['Timeout'])
            path = config['SMART CONTRACT']["Abi"]
        except:
            print('An error occured while reading config file')
        try:
            self.web3 = Web3(Web3.HTTPProvider(addr, request_kwargs={'timeout':timeout}))
            address = self.web3.toChecksumAddress(config['SMART CONTRACT']['Address'])
        except:
            print("Could not establish connection with the Blockchain test network")
        try:
            with open("Blockchain/build/contracts/Verifier.json", 'r') as fp:
                json_info = json.load(fp)
                abi = json_info['abi']
                self.account = self.web3.eth.accounts[0]
                self.web3.eth.default_account = self.account

                self.contract = self.web3.eth.contract(address=address, abi=abi)
        except:
            print('An error occured while loading abi from file')
        

    def verify(self, commitment, proof, index, value):
        return self.contract.functions.verify(commitment, proof, index, value).call()