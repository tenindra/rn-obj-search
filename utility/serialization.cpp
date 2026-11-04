// Copyright (c) 2025-2026 Tenindra Abeywickrama
// SPDX-License-Identifier: MIT

#include "serialization.h"

#include "../processing/Graph.h"
#include "../processing/Gtree.h"
#include "../processing/PackedRtree.h"
#include "../processing/ALT.h"
#include "../processing/NVDQuadtree.h"
#include "../processing/SULTree.h"
#include "../processing/HierarchicalPartition.h"

#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace serialization
{
    template <typename T>
    T getIndexFromBinaryFile(std::string idxFilePath) {
        T index;
      
        std::ifstream binaryIndexFS(idxFilePath,std::ios::in | std::ios::binary);
        if (binaryIndexFS.is_open()) {
            boost::archive::binary_iarchive ia(binaryIndexFS);
            ia >> index;
        } else {
            std::cerr << "Cannot open binary index file " << idxFilePath << std::endl;
            exit(1);
        }
        binaryIndexFS.close();
        
        return index;
    }
 
    template <typename T>
    void outputIndexToBinaryFile(T& index, std::string idxFilePath) {
        std::ofstream binaryIndexFS(idxFilePath,std::ios::out | std::ios::binary | std::ios::trunc);
        {
            boost::archive::binary_oarchive oa(binaryIndexFS);
            oa << index;
        }
        binaryIndexFS.close();
    }
    
    // Instatiate Templates
    template Graph getIndexFromBinaryFile<Graph>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Graph>(Graph& index, std::string idxFilePath);
    template Gtree getIndexFromBinaryFile<Gtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<Gtree>(Gtree& index, std::string idxFilePath);
    template OccurenceList getIndexFromBinaryFile<OccurenceList>(std::string idxFilePath);
    template void outputIndexToBinaryFile<OccurenceList>(OccurenceList& index, std::string idxFilePath);
    template ALT getIndexFromBinaryFile<ALT>(std::string idxFilePath);
    template void outputIndexToBinaryFile<ALT>(ALT& index, std::string idxFilePath);
    template ObjectList getIndexFromBinaryFile<ObjectList>(std::string idxFilePath);
    template void outputIndexToBinaryFile<ObjectList>(ObjectList& index, std::string idxFilePath);
    template PackedRtree getIndexFromBinaryFile<PackedRtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<PackedRtree>(PackedRtree& index, std::string idxFilePath);
    template NVDQuadtree getIndexFromBinaryFile<NVDQuadtree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<NVDQuadtree>(NVDQuadtree& index, std::string idxFilePath);
    template SULTree getIndexFromBinaryFile<SULTree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<SULTree>(SULTree& index, std::string idxFilePath);
    template COLTree getIndexFromBinaryFile<COLTree>(std::string idxFilePath);
    template void outputIndexToBinaryFile<COLTree>(COLTree& index, std::string idxFilePath);
    template HierarchicalPartition getIndexFromBinaryFile<HierarchicalPartition>(std::string idxFilePath);
    template void outputIndexToBinaryFile<HierarchicalPartition>(HierarchicalPartition& index, std::string idxFilePath);
}
