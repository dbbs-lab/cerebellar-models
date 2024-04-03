##################
List of strategies
##################

:class:`ConnectomeMossyGlomerulus <.connectome.to_glomerulus.ConnectomeMossyGlomerulus>`
========================================================================================
According to literature data (Billings et al., 2014; Ito, 1984; Sultan, 2001),
The algorithm selects a mf within the `60` μm along the x-axis and `20` μm along the
z-axis box surrounding each glom. This selection is random and performed with a truncated
exponential distribution. Since the placement of mf and glom is uniformed within their partition,
the convergence and divergence ratios of this connection is guaranteed.

:class:`ConnectomeGlomerulusGolgi <.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi>`
===========================================================================================

The algorithm selects all glom within the sphere (radius `50` μm) surrounding each GoC soma.
For each unique glom to connect, the tip of a basal dendrite branch from the golgi morphology is
selected. This selection is random and performed with a truncated exponential distribution based on
the distance between the tip of each branch and the glom to connect.

:class:`ConnectomeGlomerulusGranule <cerebellum.connectome.glomerulus_granule.ConnectomeGlomerulusGranule>`
===========================================================================================================

The algorithm selects `4` unique glom within the sphere (radius `40` μm) surrounding each GrC soma.
Moreover, each presynaptic glom should belong to a unique mf cluster, i.e. each should be connected
through the `ConnectomeMossyGlomerulus` strategy to a different mf. The mf cluster, the presynaptic
glom and the postsynaptic GrC dendrite are all randomly chosen. If not enough glomerulus could be
found in the `40` μm radius sphere surrounding the GrC soma, the closest glom from the remaining
cluster are selected to connect.

:class:`ConnectomeGolgiGlomerulus <cerebellum.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus>`
=====================================================================================================

The algorithm selects the closest glom (maximum `40`) are within the sphere (radius `150` μm)
surrounding each GoC soma. For each unique glom selected, the tip of an axon branch from the golgi
morphology is randomly selected. All GrC connected to the selected glom through the
`ConnectomeGolgiGlomerulus` strategy are also connected to the selected presynaptic GoC axon tip.

