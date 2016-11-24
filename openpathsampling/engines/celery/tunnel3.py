from sshtunnel import SSHTunnelForwarder

import time

def main():

    db_server = 'shark.imp.fu-berlin.de'
    redis_db_port = 6379

    redis_server = (db_server, 22)
    node_remote = ('localhost', redis_db_port)

    node_port = 6379  # can be any port. Needs to be the same as in redis worker

    redis_server_user = 'jprinz'
    keyfile = '/Users/jan-hendrikprinz/.ssh/known_hosts'

    keylines = open(keyfile).readlines()

    known_hosts_line = [l for l in keylines if l.startswith(db_server)][0]

    print known_hosts_line

    password = open('pw').read()

    server = SSHTunnelForwarder(
        db_server,
        # ssh_host_key=known_hosts_line,
        ssh_username=redis_server_user,
        ssh_password=password,
        local_bind_address=('127.0.0.1', node_port),
        remote_bind_address=('127.0.0.1', redis_db_port)
    )

    server.start()

    print(server.local_bind_port)  # show assigned local port
    # work with `SECRET SERVICE` through `server.local_bind_port`.

    time.sleep(4)

    server.stop()

if __name__ == '__main__':
    main()