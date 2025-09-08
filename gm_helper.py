#!/usr/bin/python3

import zmq
import os
from contextlib import contextmanager
import sys


@contextmanager
def zmq_client(file):
    """Creates a context manager with ZMQ client socket as resource using the IPC transport protocol, used via python with block.
    file passed as argument must be associated with a ZMQ connection openend by a server.
    """
    with zmq.Context() as ctx:
        with ctx.socket(zmq.REQ) as socket:
            addr = os.path.abspath(file)
            socket.connect("ipc://%s" % addr)
            yield socket


####################################################################
##                   script execution starts here                 ##
####################################################################

if __name__ == '__main__':
    with zmq_client(".ipc_file") as socket:
        # g16 hands over two filenames namely the input file for the external
        # and the outputfile which is the inputfile for g16.
        # Processing of the files will take place in the main script.
        infile = sys.argv[2]
        outfile = sys.argv[3]
        # send message to main script, notify to start calculation
        socket.send_string(f'{infile}|{outfile}')

        # wait for message from main script that indicates calculation is finished
        message = socket.recv_string()
