

:class:`ConnectomeMossyGlomerulus <.connectome.to_glomerulus.ConnectomeMossyGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects a mossy fiber within the `x_length` μm along the x-axis and `z_length` μm
along the z-axis box surrounding each glomerulus. This selection is random and performed with a
truncated exponential distribution. Since the placement of mossy fibers and glomerulus is uniformed
within their partition, the convergence and divergence ratios of this connection is guaranteed.

:class:`ConnectomeGlomerulusGolgi <.connectome.glomerulus_golgi.ConnectomeGlomerulusGolgi>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects all glomeruli within the sphere (radius `radius` μm) surrounding each golgi
soma. For each unique glomerulus to connect, the tip of a basal dendrite branch from the golgi
morphology is selected. This selection is random and performed with a truncated exponential
distribution based on the distance between the tip of each branch and the glomerulus to connect.

:class:`ConnectomeGlomerulusGranule <.connectome.glomerulus_granule.ConnectomeGlomerulusGranule>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects `convergence` unique glomeruli within the sphere (radius `radius` μm)
surrounding each granule cell soma. Moreover, each presynaptic glomerulus should belong to a unique
mossy fiber cluster, i.e. each should be connected through the `ConnectomeMossyGlomerulus` strategy
to a different mossy fiber. The mossy fiber cluster, the presynaptic glomerulus and the
postsynaptic granule cell dendrite are all randomly chosen. If not enough glomerulus could be
found in the `radius` μm radius sphere surrounding the granule cell soma, the closest glomerulus
from the remaining cluster are selected to connect.

:class:`ConnectomeGolgiGlomerulus <.connectome.golgi_glomerulus.ConnectomeGolgiGlomerulus>`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The algorithm selects the closest glomerulus (maximum `divergence`) are within the sphere (radius
`radius` μm) surrounding each golgi soma. For each unique glomerulus selected, the tip of an axon
branch from the golgi morphology is randomly selected. All granule cells connected to the selected
glomerulus through the `ConnectomeGolgiGlomerulus` strategy are also connected to the selected
presynaptic golgi axon tip.
