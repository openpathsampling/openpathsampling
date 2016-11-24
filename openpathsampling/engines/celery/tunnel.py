#!/usr/bin/env python

# Copyright (C) 2003-2007  Robey Pointer <robeypointer@gmail.com>
#
# This file is part of paramiko.
#
# Paramiko is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Paramiko is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Paramiko; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA.

"""
Sample script showing how to do local port forwarding over paramiko.
This script connects to the requested SSH server and sets up local port
forwarding (the openssh -L option) from a local port through a tunneled
connection to a destination reachable from the SSH server machine.
"""

import select

try:
    import SocketServer
except ImportError:
    import socketserver as SocketServer

import sys

import paramiko

SSH_PORT = 22
DEFAULT_PORT = 4000

g_verbose = True


class ForwardServer(SocketServer.ThreadingTCPServer):
    daemon_threads = True
    allow_reuse_address = True


class Handler(SocketServer.BaseRequestHandler):
    def handle(self):
        try:
            chan = self.ssh_transport.open_channel('direct-tcpip',
                                                   (self.chain_host,
                                                    self.chain_port),
                                                   self.request.getpeername())
        except Exception as e:
            verbose('Incoming request to %s:%d failed: %s' % (self.chain_host,
                                                              self.chain_port,
                                                              repr(e)))
            return
        if chan is None:
            verbose(
                'Incoming request to %s:%d was rejected by the SSH server.' %
                (self.chain_host, self.chain_port))
            return

        verbose('Connected!  Tunnel open %r -> %r -> %r' % (
        self.request.getpeername(),
        chan.getpeername(), (self.chain_host, self.chain_port)))
        while True:
            r, w, x = select.select([self.request, chan], [], [])
            if self.request in r:
                data = self.request.recv(1024)
                if len(data) == 0:
                    break
                chan.send(data)
            if chan in r:
                data = chan.recv(1024)
                if len(data) == 0:
                    break
                self.request.send(data)

        peername = self.request.getpeername()
        chan.close()
        self.request.close()
        verbose('Tunnel closed from %r' % (peername,))


def forward_tunnel(local_port, remote_host, remote_port, transport):
    # this is a little convoluted, but lets me configure things for the Handler
    # object.  (SocketServer doesn't give Handlers any way to access the outer
    # server normally.)
    class SubHander(Handler):
        chain_host = remote_host
        chain_port = remote_port
        ssh_transport = transport

    ForwardServer(('', local_port), SubHander)


def verbose(s):
    if g_verbose:
        print(s)


def main():
    # options, server, remote = parse_options()

    db_server = 'shark.imp.fu-berlin.de'
    redis_db_port = 6379

    redis_server = (db_server, 22)
    node_remote = ('localhost', redis_db_port)

    node_port = 6379  # can be any port. Needs to be the same as in redis worker

    redis_server_user = 'jprinz'
    keyfile = '/Users/jan-hendrikprinz/.ssh/known_hosts'
    look_for_keys = True

    password = None

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.WarningPolicy())

    verbose('Connecting to ssh host %s:%d ...' % (redis_server[0], redis_server[1]))
    try:
        client.connect(
            redis_server[0], redis_server[1],
            username=redis_server_user,
            key_filename=keyfile,
            look_for_keys=look_for_keys,
            password=password)

    except Exception as e:
        print('*** Failed to connect to %s:%d: %r' % (redis_server[0], redis_server[1], e))
        sys.exit(1)

    verbose('Now forwarding port %d to %s:%d ...' % (
    node_port, 'localhost', node_remote[1]))

    try:
        forward_tunnel(
            node_port, node_remote[0], node_remote[1],
                       client.get_transport())
    except KeyboardInterrupt:
        print('C-c: Port forwarding stopped.')
        sys.exit(0)


if __name__ == '__main__':
    main()