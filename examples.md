---
abstract: |
  Through a sequence of algorithmically defined transformations representing geological events such as sedimentation, folding, faulting, and magmatic intrusion, users can construct and visualize 3D geological models. This notebook documents the methodology behind applying these transformations to a spatial dataset of coordinates, with the goal of simulating the geological evolution of a hypothetical region. Each class within the notebook encapsulates a specific geological process, parameterized by real-world measurements.
---

The focus of this notebook is to implement and visualize geological transformations in a 3D space. By simulating processes such as layering, folding, and intrusion, we can gain insights into the geological history and structure of a region. The notebook details the application of a series of transformations to a set of spatial coordinates (`xyz`), representing geological features.

# Approach

The modeling approach taken here is procedural and based on the application of sequential transformations, each representing a different geological process. These transformations are applied to a 3D model of the subsurface, which is represented by a set of `xyz` coordinates. This approach allows for the step-wise construction of a complex geological history, providing a clear narrative from initial conditions to the final state.

Each transformation, whether it be a `Layer`, `Fold`, or `Dike`, is encapsulated within its class, with parameters specifying the nature of the transformation. These parameters are based on real-world geological measurements, such as _strike_, _dip_, and _rake_, and are applied to the model using mathematical operations that mimic natural processes.

The `ModelHistory` function serves as the orchestrator, taking in a sequence of transformations (`history`) and applying them to the initial `xyz` data. The events are applied in reverse order and the result is either a reverse-transformation in space or a "coloring" event that exits the algorithm. The spatial transformations are sent backwards through the geologic history, with the `xyz` vector being deformed, and an array called `data` to keep track of the geologic formation at each point in the model.

In the following sections, we will define each transformation class, apply a sequence of transformations to a model, and visualize the results.

# Implementation

The source for the basic geology calculator is in python [](./geo.py). To use the library, you can import the file:

```{code-cell} python
import geo
```

# Geologic Events

Here we define a few events.

```{code-cell} python
layer0 = geo.Layer(base=-5., width=5., value=0)
layer1 = geo.Layer(base=0., width=5., value=1)
layer2 = geo.Layer(base=5., width=5., value=2)
tilt = geo.Tilt(strike=0, dip=20)
upright_fold = geo.Fold(strike=0, dip=90)
dike = geo.Dike(strike=0, dip=60, width=3, point=[0, 0, 0], data_value=3)
```

# Layers

Lets first look at the basic X-Z cross section if we just have the layers

```{code-cell} python
xyz, X, Y, Z = geo.getModel(50)
data = geo.ModelHistory(xyz, [layer0, layer1, layer2])
geo.plotCrossSection(data, X, Y, Z)
```

# Tilting

And then a history with a tilting event?

```{code-cell} python
xyz, X, Y, Z = geo.getModel(50)
data = geo.ModelHistory(xyz, [layer0, layer1, layer2, tilt])
geo.plotCrossSection(data, X, Y, Z)
```

# Folding

```{code-cell} python
xyz, X, Y, Z = geo.getModel(50)
data = geo.ModelHistory(xyz, [layer0, layer1, layer2, tilt, upright_fold])
geo.plotCrossSection(data, X, Y, Z)
```

# With a Dike

```{code-cell} python
xyz, X, Y, Z = geo.getModel(50)
data = geo.ModelHistory(xyz, [layer0, layer1, layer2, tilt, upright_fold, dike])
geo.plotCrossSection(data, X, Y, Z)
```

The dike can come before the folding event as well, in which case it is transformed.

```{code-cell} python
xyz, X, Y, Z = geo.getModel(50)
data = geo.ModelHistory(xyz, [layer0, layer1, layer2, tilt, dike, upright_fold])
geo.plotCrossSection(data, X, Y, Z)
```
