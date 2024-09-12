import nest
from bsb import Scaffold, options, parse_configuration_file

options.verbosity = 4
cfg = parse_configuration_file("configurations/mouse/dcn-io/dcn.yaml")
network = Scaffold(cfg)

network.compile(clear=True)
