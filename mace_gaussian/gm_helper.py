#!/usr/bin/python3

import os
import sys
from contextlib import contextmanager
from pathlib import Path

import zmq


@contextmanager
def zmq_client(file):
    """Creates a context manager with ZMQ client socket as resource using IPC transport.

    The file passed as argument must be associated with a ZMQ connection opened by a server.
    Used via python with block.
    """
    with zmq.Context() as ctx, ctx.socket(zmq.REQ) as socket:
        addr = Path(file).resolve()
        socket.connect(f"ipc://{addr}")
        yield socket


####################################################################
##                   script execution starts here                 ##
####################################################################

if __name__ == "__main__":
    ipc_path = os.environ.get("MACE_IPC_PATH", ".ipc_file")
    with zmq_client(ipc_path) as socket:
        # g16 hands over two filenames namely the input file for the external
        # and the outputfile which is the inputfile for g16.
        # Processing of the files will take place in the main script.
        infile = sys.argv[2]
        outfile = sys.argv[3]
        # send message to main script, notify to start calculation
        socket.send_string(f"{infile}|{outfile}")

        # wait for message from main script that indicates calculation is finished
        message = socket.recv_string()
