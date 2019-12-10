def plot_trisurf(self, *args, **kwargs):
    """
    ============= ================================================
    Argument      Description
    ============= ================================================
    *X*, *Y*, *Z* Data values as 1D arrays
    *color*       Color of the surface patches
    *cmap*        A colormap for the surface patches.
    *norm*        An instance of Normalize to map values to colors
    *vmin*        Minimum value to map
    *vmax*        Maximum value to map
    *shade*       Whether to shade the facecolors
    ============= ================================================
    The (optional) triangulation can be specified in one of two ways;
    either::
      plot_trisurf(triangulation, ...)
    where triangulation is a :class:`~matplotlib.tri.Triangulation`
    object, or::
      plot_trisurf(X, Y, ...)
      plot_trisurf(X, Y, triangles, ...)
      plot_trisurf(X, Y, triangles=triangles, ...)
    in which case a Triangulation object will be created.  See
    :class:`~matplotlib.tri.Triangulation` for a explanation of
    these possibilities.
    The remaining arguments are::
      plot_trisurf(..., Z)
    where *Z* is the array of values to contour, one per point
    in the triangulation.
    Other arguments are passed on to
    :class:`~mpl_toolkits.mplot3d.art3d.Poly3DCollection`
    **Examples:**
    .. plot:: mpl_examples/mplot3d/trisurf3d_demo.py
    .. plot:: mpl_examples/mplot3d/trisurf3d_demo2.py
    .. versionadded:: 1.2.0
        This plotting function was added for the v1.2.0 release.
    """

    had_data = self.has_data()

    # TODO: Support custom face colours
    color = np.array(colorConverter.to_rgba(kwargs.pop('color', 'b')))

    cmap = kwargs.get('cmap', None)
    norm = kwargs.pop('norm', None)
    vmin = kwargs.pop('vmin', None)
    vmax = kwargs.pop('vmax', None)
    linewidth = kwargs.get('linewidth', None)
    shade = kwargs.pop('shade', cmap is None)
    lightsource = kwargs.pop('lightsource', None)

    tri, args, kwargs = Triangulation.get_from_args_and_kwargs(*args, **kwargs)
    z = np.asarray(args[0])

    triangles = tri.get_masked_triangles()
    xt = tri.x[triangles][...,np.newaxis]
    yt = tri.y[triangles][...,np.newaxis]
    zt = np.array(z)[triangles][...,np.newaxis]

    verts = np.concatenate((xt, yt, zt), axis=2)

    # Only need these vectors to shade if there is no cmap
    if cmap is None and shade:
        totpts = len(verts)
        v1 = np.empty((totpts, 3))
        v2 = np.empty((totpts, 3))
        # This indexes the vertex points
        which_pt = 0

    colset = []
    for i in xrange(len(verts)):
        avgzsum = verts[i,0,2] + verts[i,1,2] + verts[i,2,2]
        colset.append(avgzsum / 3.0)

        # Only need vectors to shade if no cmap
        if cmap is None and shade:
            v1[which_pt] = np.array(verts[i,0]) - np.array(verts[i,1])
            v2[which_pt] = np.array(verts[i,1]) - np.array(verts[i,2])
            which_pt += 1

    if cmap is None and shade:
        normals = np.cross(v1, v2)
    else:
        normals = []

    polyc = art3d.Poly3DCollection(verts, *args, **kwargs)

    if cmap:
        colset = np.array(colset)
        polyc.set_array(colset)
        if vmin is not None or vmax is not None:
            polyc.set_clim(vmin, vmax)
        if norm is not None:
            polyc.set_norm(norm)
    else:
        if shade:
            colset = self._shade_colors(color, normals)
        else:
            colset = color
        polyc.set_facecolors(colset)

    self.add_collection(polyc)
    self.auto_scale_xyz(tri.x, tri.y, z, had_data)

    return polyc

