# UM-Bridge
It is also possible to connect the 'SparseGridsKit.jl' to any model via the ['UM-Bridge'](https://github.com/UM-Bridge/umbridge) Julia client.

The required setup is sketched below.

```
using UMBridge
using SparseGridsKit

url = "model.url:port"

# Set up a model "modelname"
model = UMBridge.HTTPModel("modelname", url)

# Create sparse grid
n=1
k=5
mi = create_smolyak_miset(n,k)
domain = [[-1,1]]
sg = create_sparsegrid(mi,domain)

Z = get_grid_points(sg)

# Model evaluation with configuration parameters
config = Dict("Config1" => 1);

f_on_Z = [UMBridge.evaluate(model, [map(z[1])], config) for z in Z]
```