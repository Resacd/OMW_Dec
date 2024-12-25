#pragma once
#include "Node.h"
#include <fstream>
#include <string>
#include <mpi.h>
#include "Utility.h"
#include <filesystem>
#include <memory> 
#include <map>
#include <vector>
#include <string>


#ifndef _WIN32
#include <sys/stat.h>
#else
#include <windows.h>
#endif



#include "Struct.h"
#include "HaloLayer.h"
#include "Subdomain.h"

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkDelaunay3D.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkTriangleFilter.h>


#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkTriangle.h>
#include <vtkTriangle.h>

#include <vtkPolygon.h>
#include <vtkIdTypeArray.h>
#include <map>
#include <tuple>
#include <cmath>
#include <mpi.h>
#include <filesystem>
#include <sys/stat.h>
#include <string>
#include <limits>


class OutputFile
{
private:

	std::mutex m_fileMutex;
	int m_outputIdx;
	std::mutex m_dataMutex;
	Utility Utility_;

	double FindClosestLatticePoint(SubDomain& subdomain, Point& targetPoint, InterpolationAttribute attribute);
	void InterpolateAttributesToSTLVertices(SubDomain& subdomain, std::vector<Triangle>& stlTriangles, InterpolationAttribute attribute);

public:

	OutputFile() = default;
	OutputFile(int outputIdx);
	~OutputFile();


	void WriteVTK_LocalData(SubDomain &subdomain);
        void WriteVTK_Collective(SubDomain& subdomain);
	void WriteSurfaceAttribute(SubDomain &subdomai, int iter);

};

