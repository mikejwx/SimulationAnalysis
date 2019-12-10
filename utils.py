import numpy as np

##### This is an extension from Marks's cloud tracking code  ####

def dist(pos1, pos2):
    return np.sqrt((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2)


def _test_indices(i, j, diagonal=False, extended=False):
    if extended:
        # Count any cells in a 5x5 area centred on the current i, j cell as being adjacent.
        indices = []
        for ii in range(i - 2, i + 3):
            for jj in range(j - 2, j + 3):
                indices.append((ii, jj))
    else:
        # Standard, cells sharing a border are adjacent.
        indices = [(i-1, j), (i+1, j), (i, j-1), (i, j+1)]
        if diagonal:
            # Diagonal cells considered adjacent.
            indices += [(i-1, j-1), (i-1, j+1), (i+1, j-1), (i+1, j+1)]
    return indices

def cloud_edge(cld_indices, mask, wrap=True, diagonal=False):    # detect cloud edge
    cld_pos = []
    for npts in range(0, len(cld_indices[0])):
        cld_pos += [(cld_indices[0][npts], cld_indices[1][npts])]
<<<<<<< HEAD
    
=======

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
    cld_edge = []
    cld_in   = []
    for ii, jj in cld_pos:
        fault = 0
        for it, jt in _test_indices(ii, jj, diagonal):
            if not wrap:
                if it < 0 or it >= mask.shape[0] or \
                                jt < 0 or jt >= mask.shape[1]:
<<<<<<< HEAD
                    # if the point is outside the edge of the domain
=======
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
                    continue
                else:
                    it %= mask.shape[0]
                    jt %= mask.shape[1]
<<<<<<< HEAD
            
            if (it, jt) not in cld_pos:
                fault += 1
            
        if diagonal:
            if fault > 0 and fault < 8:
=======

            if (it, jt) not in cld_pos:
                fault += 1

        if diagonal:
            if fault>0 and fault <8:
>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
               cld_edge += [(ii, jj)]
            else: 
               cld_in   += [(ii, jj)]
        else:
<<<<<<< HEAD
            if fault > 0 and fault < 4:
               cld_edge += [(ii, jj)]
            else:
               cld_in   += [(ii, jj)]
        
=======
            if fault>0 and fault <4:
               cld_edge += [(ii, jj)]
            else:
               cld_in   += [(ii, jj)]

>>>>>>> 62279a4ff074ba2a906de6d583e07e7c7ce0c696
    return cld_edge, cld_in


def cloud_center(cld_indices, cld_var):    ## weighted cloud center
    cld_center_x = np.sum(cld_indices[0]*cld_var)/np.sum(cld_var)
    cld_center_y = np.sum(cld_indices[1]*cld_var)/np.sum(cld_var)

    #cld_pos = []
    #for npts in range(0, len(cld_indices[0])):
    #    cld_pos += [(cld_indices[0][npts], cld_indices[1][npts])]

    #if (cld_center_x, cld_center_y) not in cld_pos:
    #    return -1, -1
    #else:
    #    return cld_center_x, cld_center_y
    return cld_center_x, cld_center_y

def cloud_geometry_center(cld_edge, cld_in):    ## cloud geometric center
    if len(cld_in)>0:
        n_choose = 0
        dist_var = 9999999.9
    
        for npts in range(0, len(cld_in)):
            xc           = cld_in[npts][1]
            yc           = cld_in[npts][0]
            pos_c        = [yc, xc]
            in_edge_dist = []
            for n_edge in range(0, len(cld_edge)):
                x_edge       = cld_edge[n_edge][1]
                y_edge       = cld_edge[n_edge][0]
                pos_edge     = [y_edge, x_edge]
                #dis          = dist(pos_c, pos_edge)
                dis          = np.sqrt((pos_edge[0]-pos_c[0])**2 + (pos_edge[1]-pos_c[1])**2)
                in_edge_dist.append(dis)
            if np.var(in_edge_dist)<dist_var:
                dist_var = np.var(in_edge_dist)
                n_choose = npts

        gcenter_x = cld_in[n_choose][0]             
        gcenter_y = cld_in[n_choose][1]

        return gcenter_x, gcenter_y             
  
    else:
        return 0, 0

def cloud_radius(cld_center, cld_edge_pts, dx):    ## distance from cloud center to edge
    cld_radius = []
    for i in range(0, len(cld_edge_pts)):
        cld_radius.append(dist(cld_center, cld_edge_pts[i])*dx)
    
    cld_mean_radius = np.mean(cld_radius)
    cld_min_radius  = np.min(cld_radius)
    cld_max_radius  = np.min(cld_radius)

    return cld_radius, cld_mean_radius, cld_min_radius, cld_max_radius        

def label_clds(mask, diagonal=False, wrap=True, min_cells=0):
    """
    Label contiguous grid-cells with a given index - 1-max_label.
    :param np.ndarray mask: 2D mask of True/False representing (thresholded) clouds.
    :param bool diagonal: Whether to treat diagonal cells as contiguous.
    :param bool wrap: Whether to wrap on edge.
    :param int min_cells: Minimum number of grid-cells to include in a cloud.
    :return tuple(int, np.ndarray): max_label and 2D array of ints.
    """
    labels = np.zeros_like(mask, dtype=np.int32)
    max_label = 0
    acceptable_blobs = []
    for j in range(mask.shape[1]):
        for i in range(mask.shape[0]):
            if labels[i, j]:
                continue

            if mask[i, j]:
                blob_count = 1
                max_label += 1
                labels[i, j] = max_label
                outers = [(i, j)]
                while outers:
                    new_outers = []
                    for ii, jj in outers:
                        for it, jt in _test_indices(ii, jj, diagonal):
                            if not wrap:
                                if it < 0 or it >= mask.shape[0] or \
                                                jt < 0 or jt >= mask.shape[1]:
                                    continue
                            else:
                                it %= mask.shape[0]
                                jt %= mask.shape[1]

                            if not labels[it, jt] and mask[it, jt]:
                                blob_count += 1
                                new_outers.append((it, jt))
                                labels[it, jt] = max_label
                    outers = new_outers

                if blob_count >= min_cells:
                    acceptable_blobs.append(max_label)

    if min_cells > 0:
        out_blobs = np.zeros_like(labels)
        num_acceptable_blobs = 1
        for blob_index in acceptable_blobs:
            out_blobs[labels == blob_index] = num_acceptable_blobs
            num_acceptable_blobs += 1

        return num_acceptable_blobs, out_blobs
    else:
     	return max_label, labels

