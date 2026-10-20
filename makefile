core := source/Mesh.cpp source/Operators.cpp source/Solver.cpp source/SparseMatrix.cpp source/print.cpp

###### Transient heat problem
heatTransient : TransientFourier
	./TransientFourier

TransientFourier : $(core) source/mainTransFourier.cpp source/TransientFourier.cpp source/CFL.cpp
	g++ -g $(core) source/mainTransFourier.cpp source/TransientFourier.cpp source/CFL.cpp -o TransientFourier

###### Stationary heat problem
heatStationary : StationaryFourier
	./StationaryFourier

StationaryFourier : $(core) source/mainStatFourier.cpp
	g++ -g $(core) source/mainStatFourier.cpp -o StationaryFourier
