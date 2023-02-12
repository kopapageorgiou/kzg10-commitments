from brownie import accounts, Verifier
from web3 import Web3
import configparser, os
from pathlib import Path

def main():
    #admin = accounts.load('papas')
    #admin = to_address("0xd376798453D13F3EB606B9FB83833ab664C5Cbc2")
    web3 = Web3(Web3.HTTPProvider("http://127.0.0.1:8545", request_kwargs={'timeout':600}))
    #admin = web3.eth.accounts[0]
    admin = accounts[0]
    deployed = Verifier.deploy({
        "from": admin
    })
    config = configparser.ConfigParser()
    file_path = str(Path(__file__).parents[2]) + "/settings.ini"
    config.read(file_path)
    sc = config["SMART CONTRACT"]
    sc["address"] = deployed.address
    try:
        with open(file_path, 'w') as f:
            config.write(f) # Update
        f.close()
    except IOError as e:
        print(e)
    

    #print(deployed)