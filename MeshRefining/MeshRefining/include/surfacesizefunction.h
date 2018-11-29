#ifndef SURFACESIZEFUNCTION_H
#define SURFACESIZEFUNCTION_H

#include "data_io.h"

namespace sizefunction {

int opt_surface_size_function(GM3Data&gm3,VTKInput &vtk);


int input_gm3(GM3Data&gm3);

int input_vtk(VTKInput &vtk);
}
#endif // SURFACESIZEFUNCTION_H
