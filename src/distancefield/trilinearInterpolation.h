#define TRILINEAR_INTERPOLATION(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011) \
( (wx) * (wy) * (wz) *             (v111) + \
  (wx) * (wy) * (1-(wz)) *         (v110) + \
  (wx) * (1-(wy)) * (wz) *         (v101) + \
  (wx) * (1-(wy)) * (1-(wz)) *     (v100) + \
  (1-(wx)) * (wy) * (wz) *         (v011) + \
  (1-(wx)) * (wy) * (1-(wz)) *     (v010) + \
  (1-(wx)) * (1-(wy)) * (wz) *     (v001) + \
  (1-(wx)) * (1-(wy)) * (1-(wz)) * (v000))

#define GRADIENT_COMPONENT_X(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011,gridX) \
  (((wy) * (wz) *             (v111) + \
    (wy) * (1-(wz)) *         (v110) + \
    (1-(wy)) * (wz) *         (v101) + \
    (1-(wy)) * (1-(wz)) *     (v100) + \
    (-1) * (wy) * (wz) *         (v011) + \
    (-1) * (wy) * (1-(wz)) *     (v010) + \
    (-1) * (1-(wy)) * (wz) *     (v001) + \
    (-1) * (1-(wy)) * (1-(wz)) * (v000) ) / (gridX))

#define GRADIENT_COMPONENT_Y(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011,gridY) \
  (((wx) * (wz) *             (v111) + \
    (wx) * (1-(wz)) *         (v110) + \
    (wx) * (-1) * (wz) *         (v101) + \
    (wx) * (-1) * (1-(wz)) *     (v100) + \
    (1-(wx)) * (wz) *         (v011) + \
    (1-(wx)) * (1-(wz)) *     (v010) + \
    (1-(wx)) * (-1) * (wz) *     (v001) + \
    (1-(wx)) * (-1) * (1-(wz)) * (v000)) / (gridY))
  
#define GRADIENT_COMPONENT_Z(wx,wy,wz,v000,v100,v110,v010,v001,v101,v111,v011,gridZ) \
  (((wx) * (wy) *                (v111) + \
    (wx) * (wy) * (-1) *         (v110) + \
    (wx) * (1-(wy)) *            (v101) + \
    (wx) * (1-(wy)) * (-1) *     (v100) + \
    (1-(wx)) * (wy) *            (v011) + \
    (1-(wx)) * (wy) * (-1) *     (v010) + \
    (1-(wx)) * (1-(wy)) *        (v001) + \
    (1-(wx)) * (1-(wy)) * (-1) * (v000)) / (gridZ))

