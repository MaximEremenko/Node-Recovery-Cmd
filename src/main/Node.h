#ifndef NODE_H
#define NODE_H

#include "FieldElement.h"

class Node {
public:
    Node () {
    }
    
    Node (FieldElement sigma, FieldElement data[1024]) {
        m_sigma = sigma;
        for (int i = 0; i < 1024; ++i)
            m_data[i] = data[i];
    }

    FieldElement getDataSum(unsigned int* idx) {
        FieldElement res;
        for (int i = 0; i < 4; ++i) {
            res += m_data[idx[i]];
        }
        return res;
    }

    FieldElement getData(unsigned int i) {
        return m_data[i];
    }

    FieldElement* getData() {
        return m_data;
    }


    void setData(unsigned int idx, FieldElement data) {
        m_data[idx] = data;
    }
    
    const Node& operator=(const Node& node) {
        m_sigma = node.m_sigma;
        for (int i = 0; i < 1024; ++i)
            m_data[i] = node.m_data[i];
        return *this;
    }
    
private:
    FieldElement m_sigma;
    FieldElement m_data[1024];
};

#endif