import networkx as nx

class ReplicaNetwork(object):

    def __init__(self, replicas, storage=None):
        self.analysis = { } 
        self.replicas = replicas
        self.storage = storage

    def check_storage(self, storage):
        if storage != None:
            if storage != self.storage:
                self.analysis = { }
            self.storage = storage
        if self.storage == None:
            raise RuntimeError("No storage given for analysis")

    def analyze(self, storage, force=False):
        self.check_storage(storage)
        if force == False and self.analysis != { }:
            return self.analysis

        pass

    def diagram(self, storage=None, force=False):
        self.check_storage(storage)
        pass

    def flow(self, bottom, top, storage=None, force=False):
        self.check_storage(storage)
        # trips(bottom, top)
        # trips(top, bottom)
        pass

    def round_trips(self, bottom, top, storage=None, force=False):
        self.check_storage(storage)
        pass


def trace_interfaces_for_replica(replica, storage):
    pass

def trace_replicas_for_interface(interface, storage):
    pass
