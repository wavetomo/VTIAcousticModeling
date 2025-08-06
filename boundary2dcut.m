function snap = boundary2dcut(snap, nz, nx, nlayer)

snap = snap(nlayer+1:nlayer+nz, nlayer+1:nlayer+nx);

end