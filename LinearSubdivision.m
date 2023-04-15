function [newVertices,newFaces]=LinearSubdivision(vertices,faces)
newFaces = [];
newVertices = vertices;
nVertices = size(vertices,1);
nFaces = size(faces,1);
edgeVertex= zeros(nVertices, nVertices);
newIndexOfVertices = nVertices;

    for i=1 : nFaces
        [vaIndex,vbIndex,vcIndex] = deal(faces(i,1), faces(i,2),faces(i,3));
        vpIndex = addEdgeVertex( vaIndex,vbIndex) ;
        vqIndex = addEdgeVertex(vbIndex,vcIndex);
        vrIndex = addEdgeVertex(vaIndex, vcIndex);
        fourFaces = [vaIndex, vpIndex, vrIndex; vpIndex, vbIndex, vqIndex; vrIndex, vqIndex,vcIndex;vrIndex,vpIndex,vqIndex];
        newFaces = [newFaces; fourFaces];
    end
    for v1=1:nVertices-1
        for v2=v1 :nVertices
            vNIndex = edgeVertex(v1,v2);
            if (vNIndex~=0)
            newVertices( vNIndex ,: ) = 1/2*(vertices(v1, : )+vertices( v2,: ));
            end
        end
    end


    function vNIndex = addEdgeVertex(v1Index,v2Index)
    if (v1Index>v2Index) % setting: v1 <= v2
        vTmp = v1Index;
        v1Index = v2Index;v2Index = vTmp;
    end
    if (edgeVertex(v1Index, v2Index)==0)% new vertex
        newIndexOfVertices = newIndexOfVertices+1;
        edgeVertex(v1Index, v2Index) = newIndexOfVertices;
    end
    vNIndex = edgeVertex(v1Index, v2Index );
    end
end


