include $(ESMFMKFILE)

APP = main.x

SRC = ESMF_netcdf_read.f \
      user_coupler.F90  \
      user_model1.F90  \
      user_model2.F90 \
      main.F90

OBJ = ESMF_netcdf_read.o user_coupler.o user_model1.o user_model2.o main.o 

$(APP): $(OBJ)
	$(ESMF_F90COMPILER) -g -o $(APP) $(OBJ) $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKLIBS)

%.o: %.f
	$(ESMF_F90COMPILER) $(ESMF_F90COMPILEPATHS) -c -g $<
%.o: %.F90
	$(ESMF_F90COMPILER) $(ESMF_F90COMPILEPATHS) -c -g $<

install: $(APP)

clean:
	rm -f $(APP) *.o *.mod *.txt *.vtk *.o* *.e* *.nc PET*

user_model1.o: user_model1.F90 ESMF_netcdf_read.o
