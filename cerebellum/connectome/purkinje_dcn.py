import numpy as np
from bsb import ConnectionStrategy, config


@config.node
class FixedOutdegree(ConnectionStrategy):
    """
    Connect a group of postsynaptic cell types to ``outdegree`` uniformly random
    presynaptic cells from all the presynaptic cell types.
    """

    outdegree: int = config.attr(type=int, required=True)

    def connect(self, pre, post):
        out_ = self.outdegree
        rng = np.random.default_rng()
        high = sum(len(ps) for ps in post.placement)
        for ps in pre.placement:
            l = len(ps)
            pre_targets = np.full((l * out_, 3), -1)
            post_targets = np.full((l * out_, 3), -1)
            ptr = 0
            for i in range(l):
                pre_targets[ptr : ptr + out_, 0] = i
                post_targets[ptr : ptr + out_, 0] = rng.choice(high, out_, replace=False)
                ptr += out_
            lowmux = 0
            for post_ps in post.placement:
                highmux = lowmux + len(post_ps)
                demux_idx = (post_targets[:, 0] >= lowmux) & (post_targets[:, 0] < highmux)
                demuxed = post_targets[demux_idx]
                demuxed[:, 0] -= lowmux
                self.connect_cells(ps, post_ps, pre_targets[demux_idx], demuxed)
                lowmux = highmux
