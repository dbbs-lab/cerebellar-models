from imp import get_tag
import unittest
from bsb.unittest import NumpyTestCase
from bsb.config import from_json
from bsb.core import Scaffold
from bsb.storage import Storage
import numpy as np

# RandomStorageFixture, NumpyTestCase,
class TestGlomerulusGolgi(unittest.TestCase, NumpyTestCase):
    def setUp(self):
        self.storage = Storage(engine="hdf5", root="test")
        # self.storage = self.random_storage(engine="hdf5")
        self.cfg = from_json("cerebellum.json")
        self.cfg.network.x = 100
        self.cfg.network.y = 100
        self.cfg.network.z = 100
        self.network = Scaffold(self.cfg, self.storage)
        # self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","purkinje_layer_placement","molecular_layer_placement","mossy_fibers_to_glomerulus","glomerulus_to_golgi","glomerulus_to_granule","golgi_to_glomerulus","golgi_to_golgi","parallel_fiber_to_purkinje"])
        # self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","mossy_fibers_to_glomerulus","glomerulus_to_granule"])

    def test_convergence(self):
        cs = self.network.get_connectivity_set(tag="glomerulus_to_golgi_granule")
        # cs = self.network.get_connectivity_set(tag="golgi_cell_to_golgi_cell")

        # MOSSY-GLOMS

        convergence_list = []
        for lchunk, itr in cs.nested_iter_connections("inc"):
            # print(lchunk)
            for gchunk, (golgi_locs, glomerulus_locs) in itr:
                unique_golgi = np.unique(golgi_locs, axis=0)
                # print(unique_golgi)
                first_global = True
                n_gloms = 0
                n_golgi = 0
                for golgi in unique_golgi:
                    connected_gloms_ids = np.where(golgi_locs[:, 0] == golgi[0])
                    # print(connected_gloms_ids)
                    if first_global == True:
                        n_golgi += len(unique_golgi)
                        first_global = False
                    n_gloms += len(np.unique(connected_gloms_ids))
                    # print(n_gloms)
                # print(n_golgi)
                convergence_list.append(n_gloms / n_golgi)
                # print(len(locs))
        # print(convergence_list)
        variance = np.var(convergence_list)
        sum = 0
        for elem in convergence_list:
            sum += elem
        print(
            "Mean:",
            round(np.mean(convergence_list), 2),
            " STD:",
            round(np.sqrt(variance), 1),
        )

    def test_divergence(self):
        # cs = self.network.get_connectivity_set(tag="glomerulus_to_golgi_cell")
        cs = self.network.get_connectivity_set(tag="mossy_fibers_to_glomerulus")

        # MOSSY-GLOMS

        count_connected_glom = 0
        count_connected_mossy = 0
        list_glom = []
        list_mossy = []
        list_div = []
        for dir_, itr in cs.nested_iter_connections(direction="out"):
            first = True
            local_gloms = 0
            local_mossy = 0
            for lchunk, itr in itr:
                if first == True:
                    count_connected_mossy += len(np.unique(itr[1], axis=0))
                    print(count_connected_mossy)
                    local_mossy += len(np.unique(itr[1], axis=0))
                    first = False
                count_connected_glom += len(itr[0])
                local_gloms += len(itr[0])
            list_div.append(local_gloms / local_mossy)

        var = np.var(list_div)
        mean = np.mean(list_div)
        print("\n")
        print("\n")
        print("----------------------------------------")
        print("Connectivity: mossy_fibers_to_glomerulus")
        print("Divergence mean:", round(mean, 1), "STD:", np.sqrt(var))

        print("Connected glomeruli:", count_connected_glom)
        ps = self.network.get_placement_set("glomerulus")
        print("Total glomeruli: ", len(ps))
        ps = self.network.get_placement_set("mossy_fibers")
        print("Connected mossy:", count_connected_mossy)
        print("Total mossy: ", len(ps))
        print("----------------------------------------")


class TestMossyGlomerulus(unittest.TestCase, NumpyTestCase):
    def setUp(self):
        self.storage = Storage(engine="hdf5", root="test_str")
        # self.storage = self.random_storage(engine="hdf5")
        self.cfg = from_json("cerebellum.json")
        self.cfg.network.x = 100
        self.cfg.network.y = 100
        self.cfg.network.z = 100
        self.network = Scaffold(self.cfg, self.storage)
        # self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","purkinje_layer_placement","molecular_layer_placement","mossy_fibers_to_glomerulus","glomerulus_to_golgi","glomerulus_to_granule","golgi_to_glomerulus","golgi_to_golgi","parallel_fiber_to_purkinje"])
        self.network.compile(
            only=[
                "granular_layer_placement",
                "granular_layer_innervation",
                "mossy_fibers_to_glomerulus",
                "glomerulus_to_golgi",
                "glomerulus_to_granule",
                "golgi_to_glomerulus",
                "golgi_to_golgi",
            ]
        )
        # self.network.compile(only=["granular_layer_placement", "granular_layer_innervation","mossy_fibers_to_glomerulus"])

    def test_convergence(self):
        cs = self.network.get_connectivity_set(tag="glomerulus_to_golgi_granule")
        # cs = self.network.get_connectivity_set(tag="golgi_cell_to_golgi_cell")

        # MOSSY-GLOMS

        convergence_list = []
        for lchunk, itr in cs.nested_iter_connections("inc"):
            # print(lchunk)
            for gchunk, (golgi_locs, glomerulus_locs) in itr:
                unique_golgi = np.unique(golgi_locs, axis=0)
                # print(unique_golgi)
                first_global = True
                n_gloms = 0
                n_golgi = 0
                for golgi in unique_golgi:
                    connected_gloms_ids = np.where(golgi_locs[:, 0] == golgi[0])
                    # print(connected_gloms_ids)
                    if first_global == True:
                        n_golgi += len(unique_golgi)
                        first_global = False
                    n_gloms += len(np.unique(connected_gloms_ids))
                    # print(n_gloms)
                # print(n_golgi)
                convergence_list.append(n_gloms / n_golgi)
                # print(len(locs))
        # print(convergence_list)
        variance = np.var(convergence_list)
        sum = 0
        for elem in convergence_list:
            sum += elem
        print(
            "Mean:",
            round(np.mean(convergence_list), 2),
            " STD:",
            round(np.sqrt(variance), 1),
        )

    def test_divergence(self):
        # cs = self.network.get_connectivity_set(tag="glomerulus_to_golgi_cell")
        cs = self.network.get_connectivity_set(tag="mossy_fibers_to_glomerulus")

        # MOSSY-GLOMS

        count_connected_glom = 0
        count_connected_mossy = 0
        list_glom = []
        list_mossy = []
        list_div = []
        for dir_, itr in cs.nested_iter_connections(direction="out"):
            first = True
            local_gloms = 0
            local_mossy = 0
            for lchunk, itr in itr:
                if first == True:
                    count_connected_mossy += len(np.unique(itr[0], axis=0))
                    print(count_connected_mossy)
                    local_mossy += len(np.unique(itr[0], axis=0))
                    first = False
                count_connected_glom += len(itr[1])
                local_gloms += len(itr[1])
            list_div.append(local_gloms / local_mossy)

        var = np.var(list_div)
        mean = np.mean(list_div)
        print("\n")
        print("\n")
        print("----------------------------------------")
        print("Connectivity: mossy_fibers_to_glomerulus")
        print("Divergence mean:", round(mean, 1), "STD:", np.sqrt(var))

        print("Connected glomeruli:", count_connected_glom)
        ps = self.network.get_placement_set("glomerulus")
        print("Total glomeruli: ", len(ps))
        ps = self.network.get_placement_set("mossy_fibers")
        print("Connected mossy:", count_connected_mossy)
        print("Total mossy: ", len(ps))
        print("----------------------------------------")
