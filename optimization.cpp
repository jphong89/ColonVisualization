#include "optimization.h"

vtkSmartPointer<vtkPolyData> Optimize(vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp, bool *Is_Fixed)
{
    vtkSmartPointer<vtkPolyData> OptimizedSurface = vtkSmartPointer<vtkPolyData>::New();
    OptimizedSurface->DeepCopy(SurfaceLineUp);
    // maintain a list to record the correspondence between original ids and the variable ids
    vtkSmartPointer<vtkIdList> Ids = vtkSmartPointer<vtkIdList>::New();
    int* InvertIds = (int*)malloc(t_colon->GetNumberOfPoints()*sizeof(int));
    for(int i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        InvertIds[i] = -1;
    }
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for(vtkIdType i=0; i < t_colon->GetNumberOfPoints(); i++)
    {
        if(!Is_Fixed[i])
        {
            double p[3];
            t_colon->GetPoint(i, p);
            Ids->InsertNextId(i);
            points->InsertNextPoint(p);
            InvertIds[i] = Ids->GetNumberOfIds() - 1;
        }
    }


    int m = Ids->GetNumberOfIds() * 3;
    std::cout<<"total number of variables = "<<m<<endl;
    std::map<vtkIdType, double> coefficientMap;
    coefficientMap.clear();

    double* b_c = (double *)malloc(m*sizeof(double));
    memset(b_c, 0, sizeof(double)*m);

    constructAandb(coefficientMap, b_c, t_colon, SurfaceLineUp, Is_Fixed, Ids, InvertIds);

    // transform coefficientMap into a data structure that eigen recognizes
    std::vector<T> coefficients;
    coefficients.clear();
    for(std::map<vtkIdType, double>::iterator it = coefficientMap.begin(); it!=coefficientMap.end(); it++)
    {
        int idx1 = (vtkIdType)it->first % (vtkIdType)m;
        int idx2 = floor((double(it->first) + 0.01) / m);

        coefficients.push_back(T(idx1, idx2, it->second));
        //if(idx1 == 51242 || idx2 == 51242)
        std::cout<<it->first<<" "<<idx1<<" "<<idx2<<" "<<it->second<<"- "<<m<<endl;
    }
    SpMat A(m,m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    Eigen::VectorXd b(m);
    b.setZero();
    for(int k = 0; k < m; k++)
    {
        b(k) = b_c[k];
    }

    Eigen::SimplicialCholesky<SpMat> chol(A);
    Eigen::VectorXd solved_x = chol.solve(b);

    for(int k = 0; k < m; k++)
    {
        double x;
        x = solved_x(k);
        std::cout<<x<<endl;
    }

    free(b_c);
    free(InvertIds);
    return OptimizedSurface;
}

void updateA(int m, int i1, int i2, double weight,std::map<vtkIdType, double> &coefficientMap)
{
    m = (vtkIdType)m;
    i1 = (vtkIdType)i1;
    i2 = (vtkIdType)i2;
    //if(i1 == 17080 && i2 == 17080)
    //std::cout<<"updatea "<<i1<<"-"<<i2<<endl;
    std::map<vtkIdType, double>::iterator iter;
    vtkIdType k = (i1*(vtkIdType)3) + (i2*(vtkIdType)3) * m;
    //std::cout<<k<<" ";

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*(vtkIdType)3+(vtkIdType)1) + (i2*(vtkIdType)3+(vtkIdType)1) * m;
    //std::cout<<k<<" ";

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }

    k = (i1*(vtkIdType)3+(vtkIdType)2) + (i2*(vtkIdType)3+(vtkIdType)2) * m;
    //if(i1 == 17080 && i2 == 17080)
    //std::cout<<k<<"- "<<i1*(vtkIdType)3+(vtkIdType)2<<" "<<i2*(vtkIdType)+(vtkIdType)2<<endl;

    iter = coefficientMap.find(k);
    if (iter == coefficientMap.end()){
        coefficientMap.insert(std::pair<vtkIdType,double>(k,weight));
    }else{
        coefficientMap[k] = iter->second + weight;
    }
}
vtkSmartPointer<vtkIdList> GetFacetsOfEdge(vtkSmartPointer<vtkPolyData> mesh, int idx1, int idx2)
{
    vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
    vtkSmartPointer<vtkIdList> list1 = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(idx1, list1);
    vtkSmartPointer<vtkIdList> list2 = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(idx2, list2);

    for(int j = 0; j < list1->GetNumberOfIds(); j++)
    {
        for(int k = 0; k < list2->GetNumberOfIds(); k++)
        {
            if(list1->GetId(j) == list2->GetId(k))
            {
                list->InsertNextId(list1->GetId(j));
            }
        }
    }
    return list;
}

double ComputeStretchWeight(vtkSmartPointer<vtkPolyData> t_colon, int idx1, int idx2)
{
    vtkSmartPointer<vtkIdList> list = vtkSmartPointer<vtkIdList>::New();
    list = GetFacetsOfEdge(t_colon, idx1, idx2);

    double area = 0;
    if(list->GetNumberOfIds() == 1 || list->GetNumberOfIds() == 2)
    {
        for(int i=0; i < list->GetNumberOfIds(); i++)
        {
            area += ComputeFacetArea(t_colon, list->GetId(i));
        }
    }
    else
    {
        std::cerr<<"error, list number not correct: "<<list->GetNumberOfIds()<<endl;
        exit(1);
    }
    double p1[3], p2[3];
    t_colon->GetPoint(idx1, p1);
    t_colon->GetPoint(idx2, p2);
    double l = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    //std::cout<<list->GetNumberOfIds()<<" - "<<area<<" "<<l<<endl;
    return area / l / l;
}

void constructAandb(std::map<vtkIdType, double> &coefficientMap, double *b,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp,
                bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids, int* InvertIds)
{
    int m = Ids->GetNumberOfIds() * 3;

    // stretch
    // extract the edges
    bool * marked = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(marked, 0, t_colon->GetNumberOfPoints()*sizeof(bool));
    bool * boundary = (bool*)malloc(t_colon->GetNumberOfPoints()*sizeof(bool));
    memset(boundary, 0, t_colon->GetNumberOfPoints()*sizeof(bool));

    vtkSmartPointer<vtkIntArray> lines = vtkSmartPointer<vtkIntArray>::New();
    lines->SetNumberOfComponents(2);
    for(int i=0; i<t_colon->GetNumberOfPoints(); i++)
    {
        vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
        t_colon->GetPointCells(i, cellids);
        for(int j=0; j<cellids->GetNumberOfIds(); j++)
        {
            vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
            t_colon->GetCellPoints(cellids->GetId(j), ptids);
            assert(ptids->GetNumberOfIds() <= 3);
            for(int k = 0; k < ptids->GetNumberOfIds(); k++)
            {
                int ptid = ptids->GetId(k);
                if(ptid == i || marked[ptid])
                    continue;
                int tuple[2];
                tuple[0] = i;
                tuple[1] = ptid;
                lines->InsertNextTypedTuple(tuple);
            }
        }
        marked[i] = true;
    }
    for(int i=0; i < lines->GetNumberOfTuples(); i++)
    {
        int tuple[2];
        lines->GetTypedTuple(i, tuple);
        int idx1 = tuple[0];
        int idx2 = tuple[1];
        //std::cout<<InvertIds[idx1]<<" "<<InvertIds[idx2]<<endl;
        if(GetFacetsOfEdge(t_colon, idx1, idx2)->GetNumberOfIds() == 1)
        {
            boundary[idx1] = true;
            boundary[idx2] = true;
        }

        if(!Is_Fixed[idx1] && !Is_Fixed[idx2]) // if the two points are all flexible, we only need to update A
        {
            int vidx1 = InvertIds[idx1];
            int vidx2 = InvertIds[idx2];
            //std::cout<<m<<" "<<vidx1<<" "<<vidx2<<endl;
            assert(vidx1 >= -0.5 && vidx2 >= -0.5);
            //std::cout<<vidx1<<","<<vidx2<<" - "<<idx1<<","<<idx2<<endl;
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            //std::cout<<"stretch weight = "<<stretchWeight<<endl;
            double weight = REGWEIGHT * 2 * stretchWeight * 2;
            updateA(m, vidx1, vidx1, weight, coefficientMap);
            updateA(m, vidx2, vidx2, weight, coefficientMap);
            updateA(m, vidx1, vidx2, -weight, coefficientMap);
            updateA(m, vidx2, vidx1, -weight, coefficientMap);
        }
        else if(!Is_Fixed[idx1]) // if one point is fixed, we need to update A for the flexible point and add a linear term in b
        {
            int vidx1 = InvertIds[idx1];
            //std::cout<<m<<" "<<vidx1<<endl;
            assert(vidx1 > -0.5);
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            double weight = REGWEIGHT * 2 * stretchWeight * 2;
            updateA(m, vidx1, vidx1, weight, coefficientMap);
            double v[3], pold[3], pnew[3];
            t_colon->GetPoint(idx2, pold);
            SurfaceLineUp->GetPoint(idx2, pnew);
            vtkMath::Subtract(pnew, pold, v);
            b[vidx1*3] +=   -2 * v[0] * weight;
            b[vidx1*3+1] += -2 * v[1] * weight;
            b[vidx1*3+2] += -2 * v[2] * weight;
        }
        else if(!Is_Fixed[idx2])
        {
            int vidx2 = InvertIds[idx2];
            //std::cout<<m<<" "<<vidx2<<endl;
            assert(vidx2 > -0.5);
            double stretchWeight = ComputeStretchWeight(t_colon, idx1, idx2);
            double weight = REGWEIGHT * 2 * stretchWeight * 2;
            updateA(m, vidx2, vidx2, weight, coefficientMap);
            double v[3], pold[3], pnew[3];
            t_colon->GetPoint(idx1, pold);
            SurfaceLineUp->GetPoint(idx1, pnew);
            vtkMath::Subtract(pnew, pold, v);
            b[vidx2*3] +=   -2 * v[0] * weight;
            b[vidx2*3+1] += -2 * v[1] * weight;
            b[vidx2*3+2] += -2 * v[2] * weight;
        }

    }

    // bend
    // compute surround areas
    vtkSmartPointer<vtkDoubleArray> SurroundAreas = vtkSmartPointer<vtkDoubleArray>::New();
    for(int i = 0; i < t_colon->GetNumberOfPoints(); i++)
    {
        double surroundarea = 0;
        vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
        t_colon->GetPointCells(i, cellids);
        for(int j = 0; j < cellids->GetNumberOfIds(); j++)
        {
            surroundarea += ComputeFacetArea(t_colon, cellids->GetId(j));
        }
        SurroundAreas->InsertNextValue(surroundarea);
    }

    for(int i = 0; i < Ids->GetNumberOfIds(); i++)
    {
        int idx1 = Ids->GetId(i);
        assert(!Is_Fixed[idx1]);
        if(boundary[idx1]) // don't process boundary edges
            continue;

        //std::cout<<i<<" - "<<idx1<<endl;
        double weight = 2* BENDWEIGHT / (SurroundAreas->GetValue(idx1) * SurroundAreas->GetValue(idx1));
        //double weight = 2* BENDWEIGHT / (SurroundAreas->GetValue(idx1)/2);
        double center_weight = 0;
        vtkSmartPointer<vtkIdList> neighbours = vtkSmartPointer<vtkIdList>::New();
        neighbours = GetConnectedVertices(t_colon, idx1);
        vtkSmartPointer<vtkDoubleArray> neighbor_weights = vtkSmartPointer<vtkDoubleArray>::New();
        for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
        {
            int idx2 = neighbours->GetId(j);
            double neighbor_weight = 0;
            // compute the weight[j] for the j-th neighbour
            vtkSmartPointer<vtkIdList> cellids = vtkSmartPointer<vtkIdList>::New();
            cellids = GetFacetsOfEdge(t_colon, idx1, idx2);
            assert(cellids->GetNumberOfIds() == 2); // we have already excluded the boundary edges, so this edge must be inside the surface

            for(int k = 0; k < 2; k++)
            {
                int idx12;
                vtkSmartPointer<vtkIdList> ptids1 = vtkSmartPointer<vtkIdList>::New();
                t_colon->GetCellPoints(cellids->GetId(0), ptids1);
                for(int k = 0; k < ptids1->GetNumberOfIds(); k++)
                {
                    if(ptids1->GetId(k)!=idx1 && ptids1->GetId(k)!=idx2)
                    {
                        idx12 = ptids1->GetId(k);
                        break;
                    }
                }
                double p1[3], p2[3], p12[3];
                t_colon->GetPoint(idx1, p1);
                t_colon->GetPoint(idx2, p2);
                t_colon->GetPoint(idx12, p12);
                double v1[3], v2[3];
                vtkMath::Subtract(p1, p12, v1);
                vtkMath::Subtract(p2, p12, v2);
                double angle = vtkMath::AngleBetweenVectors(v1, v2);

                neighbor_weight += cos(angle)/sin(angle);
            }
            neighbor_weight = (neighbor_weight > 0.1)? neighbor_weight : 0.1;
            //std::cout<<idx2<<" "<<neighbor_weight<<" | ";
            center_weight += neighbor_weight;
            neighbor_weights->InsertNextValue(neighbor_weight);
        }
        //std::cout<<endl;
        int vidx1 = InvertIds[idx1];
        updateA(m, vidx1, vidx1, weight* center_weight * center_weight, coefficientMap);

        for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
        {
            int idx2 = neighbours->GetId(j);
            int vidx2 = InvertIds[idx2];
            if(!Is_Fixed[idx2])
            {
                updateA(m, vidx1, vidx2, -weight*center_weight*neighbor_weights->GetValue(j), coefficientMap);
                updateA(m, vidx2, vidx1, -weight*center_weight*neighbor_weights->GetValue(j), coefficientMap);
            }
            else
            {
                double v[3], pold[3], pnew[3];
                t_colon->GetPoint(idx2, pold);
                SurfaceLineUp->GetPoint(idx2, pnew);
                vtkMath::Subtract(pold, pnew, v);

                // update the b of idx1
                b[vidx1*3] += -2 * center_weight * neighbor_weights->GetValue(j) * v[0];
                b[vidx1*3+1] += -2 * center_weight * neighbor_weights->GetValue(j) * v[1];
                b[vidx1*3+2] += -2 * center_weight * neighbor_weights->GetValue(j) * v[2];

                for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
                {
                    int currentidx = neighbours->GetId(k);
                    if(!Is_Fixed[currentidx])
                    {
                        int currentvidx = InvertIds[currentidx];
                        assert(currentvidx > -0.5);
                        b[currentvidx*3] = 2 * neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[0];
                        b[currentvidx*3+1] = 2 * neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[1];
                        b[currentvidx*3+2] = 2 * neighbor_weights->GetValue(k) * neighbor_weights->GetValue(j) * v[2];
                    }
                }
            }
        }
        for(int j = 0; j < neighbours->GetNumberOfIds(); j++)
        {
            for(int k = 0; k < neighbours->GetNumberOfIds(); k++)
            {
                int idx1 = neighbours->GetId(j);
                int idx2 = neighbours->GetId(k);
                int vidx1 = InvertIds[idx1];
                int vidx2 = InvertIds[idx2];
                if(!Is_Fixed[idx1] && !Is_Fixed[idx2])
                    updateA(m, vidx1, vidx2, weight*neighbor_weights->GetValue(j)*neighbor_weights->GetValue(k), coefficientMap);
            }
        }
    }

    free(marked);
    free(boundary);
}

void constructb(int idx, double *b,
                vtkSmartPointer<vtkPolyData> t_colon, vtkSmartPointer<vtkPolyData> SurfaceLineUp,
                bool *Is_Fixed, vtkSmartPointer<vtkIdList> Ids)
{

}

vtkSmartPointer<vtkIdList> GetConnectedVertices(vtkSmartPointer<vtkPolyData> mesh, vtkIdType id)
{
    vtkSmartPointer<vtkIdList> connectedVertices =
            vtkSmartPointer<vtkIdList>::New();

    //get all cells that vertex 'id' is a part of
    vtkSmartPointer<vtkIdList> cellIdList =
            vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(id, cellIdList);

    for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
    {
        vtkSmartPointer<vtkIdList> pointIdList =
                vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);

        for(vtkIdType j = 0; j < pointIdList->GetNumberOfIds(); j++)
        {
            if(pointIdList->GetId(j) != id)
            {
                bool marked = false;
                for(int k = 0; k < connectedVertices->GetNumberOfIds(); k++)
                {
                    if(connectedVertices->GetId(k) == pointIdList->GetId(j))
                    {
                        marked = true;
                        break;
                    }
                }
                if(!marked)
                    connectedVertices->InsertNextId(pointIdList->GetId(j));
            }
        }
    }
    return connectedVertices;
}
double ComputeFacetArea(vtkSmartPointer<vtkPolyData> t_colon, vtkIdType cellid)
{
    vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
    t_colon->GetCellPoints(cellid, ptids);
    assert(ptids->GetNumberOfIds() == 3);
    double p0[3], p1[3], p2[3];
    t_colon->GetPoint(ptids->GetId(0), p0);
    t_colon->GetPoint(ptids->GetId(1), p1);
    t_colon->GetPoint(ptids->GetId(2), p2);

    double a = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
    double b = sqrt(vtkMath::Distance2BetweenPoints(p0, p2));
    double c = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    double s = (a+b+c)/2;
    return sqrt(s*(s-a)*(s-b)*(s-c));
}
