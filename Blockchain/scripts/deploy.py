from ape import accounts, project, networks
import configparser, os
from pathlib import Path

def main():
    admin = accounts.load('ganache-local')
    deployed = project.Verifier.deploy(sender=admin)

    config = configparser.ConfigParser()
    file_path = str(Path(__file__).parents[2]) + "/settings.ini"
    if os.path.exists(file_path):
        config.read(file_path)
        sc = config["SMART CONTRACT"]
        sc["address"] = str(deployed)

    else:
        parts = networks.active_provider.provider_settings['uri'].split(':')
        port = parts[2]
        server_address = parts[0] + ':' + parts[1]
        config['SERVER'] = {'address': server_address, 'port': port, 'timeout': 600}
        sol_path = str(project.Verifier.source_path)
        parent_path = sol_path.split('contracts')[0]
        file_name = sol_path.split('contracts')[1].split('.')[0]
        abi_path = parent_path + '.build' + file_name + '.json'
        config['SMART CONTRACT'] = {'address': str(deployed), 'abi': abi_path}
    try:
        with open(file_path, 'w') as f:
            config.write(f) # Update
        f.close()
    except IOError as e:
        print(e)