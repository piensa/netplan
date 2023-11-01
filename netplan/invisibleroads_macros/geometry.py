def transform_geometries(geometries, f):
    transformed_geometries = []

    def transform_coords(coords):
        return [f(xyz) for xyz in coords]

    for x in geometries:
        GeometryClass = x.__class__
        if hasattr(x, 'geoms'):
            y = GeometryClass(transform_geometries(x.geoms, f))
        elif hasattr(x, 'exterior'):
            y = GeometryClass(transform_coords(x.exterior.coords), [
                transform_coords(i.coords) for i in x.interiors])
        else:
            y = GeometryClass(transform_coords(x.coords))
        transformed_geometries.append(y)
    return transformed_geometries


def flip_xy(xyz):
    'Flip x and y coordinates whether or not there is a z-coordinate'
    xyz = list(xyz)  # Preserve original
    xyz[0], xyz[1] = xyz[1], xyz[0]
    return tuple(xyz)


def drop_z(xyz):
    return tuple(xyz[:2])
