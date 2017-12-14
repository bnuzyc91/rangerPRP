# rangerPRP

//init
void findBestSplitValueLR::init(Data* data,  std::vector<std::vector<size_t>>& sampleIDs, size_t nodeID ,size_t varID, std::vector<double>* unique_timepoints,
                                size_t status_varID,std::vector<size_t>* response_timepointIDs,size_t min_node_size,size_t* num_deaths, size_t* num_samples_at_risk)
{
    
//         this->data=0;
//        this->sampleIDs=0;
//         this->nodeID= 0;
//        this->varID=0;
//        this->unique_timepoints=0;
//        this->status_varID=0;
//        this->response_timepointIDs=0;
//        this->min_node_size=0;
//   num_timepoints = unique_timepoints->size();
//        this->num_deaths =0 ;
//        this->num_samples_at_risk = 0;
    
     this->data=data;
//     this->sampleIDs=sampleIDs;
     this->nodeID= nodeID;
    this->varID=varID;
    this->unique_timepoints=unique_timepoints;
    this->status_varID=status_varID;
    this->response_timepointIDs=response_timepointIDs;
    this->min_node_size=min_node_size;
//    num_timepoints = unique_timepoints->size();
    this->num_deaths =num_deaths ;
    this->num_samples_at_risk = num_samples_at_risk;
    double value =data->get(nodeID, varID);
    size_t i=0;
    std:: cout << "inside findBestSplitValueLR value is " << value << "\n"<< std::endl;
//    std:: cout << "sampleIDs is " << sampleIDs[nodeID].size() << "\n"<< std::endl;
    std:: cout << "unique_timepoints is " << unique_timepoints[i].size() << "\n"<< std::endl;
//    std:: cout << "response_timepointIDs is " << response_timepointIDs[i].size() << "\n"<< std::endl;
//    std:: cout << "status_varID is "<< status_varID << "\n"<< std::endl;
//    std:: cout << "num_deaths is " << num_deaths[i] << "\n"<< std::endl;
//    std:: cout << "num_samples_at_risk is " << num_samples_at_risk[i] << "\n"<< std::endl;
//    std:: cout << "min_node_size is "<< min_node_size << "\n"<< std::endl;
//
//     this->varID=varID;
//    std::cout << typeid(varID).name() << std::endl;
//    std:: cout << "inside findBestSplitValueLR varID is " << varID << "\n"<< std::endl;
}
