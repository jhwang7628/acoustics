#ifndef MODAL_MATERIAL_LIST_H 
#define MODAL_MATERIAL_LIST_H 
#include <vector>

//##############################################################################
//##############################################################################
typedef std::vector<std::shared_ptr<ModalMaterial> > ModalMaterialList; 
//struct ModalMaterialList
//{
//    std::vector<std::shared_ptr<ModalMaterial> > modalMaterials; 
//    std::shared_ptr<ModalMaterial> operator[](const int &ind){return modalMaterials.at(ind);}
//    std::shared_ptr<ModalMaterial> operator()(const int &ind){return modalMaterials.at(ind);}
//    void push_back(std::shared_ptr<ModalMaterial> &material){modalMaterials.push_back(material);}
//};
#endif
