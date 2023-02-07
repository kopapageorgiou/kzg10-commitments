from brownie import accounts, Verifier
from brownie.convert import to_address
from web3 import Web3

def main():
    #admin = accounts.load('papas')
    #admin = to_address("0xd376798453D13F3EB606B9FB83833ab664C5Cbc2")
    admin = accounts[0]
    deployed = Verifier.deploy({
        "from": admin
    })
    

    #print(deployed)