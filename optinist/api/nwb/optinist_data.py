from pynwb import get_class, load_namespaces
from pynwb.spec import NWBDatasetSpec, NWBGroupSpec, NWBNamespaceBuilder

name = "optinist"
ns_path = f"{name}.namespace.yaml"
ext_source = f"{name}.extensions.yaml"

# Now we define the data structures. We use `NWBDataInterface` as the base type,
# which is the most primitive type you are likely to use as a base. The name of the
# class is `CorticalSurface`, and it requires two matrices, `vertices` and
# `faces`.

postprocess = NWBGroupSpec(
    doc="postprocess",
    datasets=[
        NWBDatasetSpec(
            doc="data",
            shape=[
                (None,),
                (None, None),
                (None, None, None),
                (None, None, None, None),
            ],
            name="data",
            dtype="float",
        )
    ],
    neurodata_type_def="PostProcess",
    neurodata_type_inc="NWBDataInterface",
)

# Now we set up the builder and add this object

ns_builder = NWBNamespaceBuilder(f"{name} extensions", name, version="0.1.0")
ns_builder.add_spec(ext_source, postprocess)
ns_builder.export(ns_path)

load_namespaces(ns_path)
PostProcess = get_class("PostProcess", name)
