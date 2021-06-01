#ifndef FIELDELEMENT_H
#define FIELDELEMENT_H

#include "FieldDefs.h"

class FieldElement {
public:
    FieldElement()
    {
        m_iElement = 0;
    }
        
    FieldElement(fe_type iElement)
    {
        m_iElement = iElement;
    }
    
    FieldElement(const FieldElement& element)
    {
        m_iElement = element.m_iElement;
    }
    
    ~FieldElement()
    {
    }
    
    fe_type getElement()
    {
        return m_iElement;
    }
    
    inline const FieldElement& operator= (const FieldElement& element)
    {
        m_iElement = element.m_iElement;
        return *this;
    }
    
    inline const FieldElement& operator= (const fe_type iElement)
    {
        m_iElement = iElement;
        return *this;
    }
    
    inline bool operator == (const FieldElement& element) const
    {
        return (m_iElement == element.m_iElement);
    }
    
    inline FieldElement operator+ (const FieldElement& element) const
    {
        FieldElement result(m_iElement^element.m_iElement);
        return result;
    }
    
    inline FieldElement operator+ (const fe_type iElement) const
    {
        FieldElement result(m_iElement^iElement);
        return result;
    }
    
    inline const FieldElement& operator+= (const FieldElement& element)
    {
        m_iElement ^= element.m_iElement;
        return *this;
    }
    
    inline const FieldElement& operator+= (const fe_type iElement)
    {
        m_iElement ^= iElement;
        return *this;
    }
    
    inline FieldElement operator* (const FieldElement& element) const
    {
        FieldElement result;

        if (m_iElement && element.m_iElement)
        {
            fe_type a_deg, b_deg, res = 0;
            a_deg = EL_TO_DEG[m_iElement - 1];
            b_deg = EL_TO_DEG[element.m_iElement - 1];
            res = (a_deg + b_deg);
            result.m_iElement = DEG_TO_EL[res - (res >= 4095)*4095];
        }

        return result;
    }
    
    inline FieldElement operator* (const fe_type iElement) const
    {
        FieldElement result;

        if (m_iElement && iElement)
        {
            fe_type a_deg, b_deg, res = 0;
            a_deg = EL_TO_DEG[m_iElement - 1];
            b_deg = EL_TO_DEG[iElement - 1];
            res = (a_deg + b_deg);
            result.m_iElement = DEG_TO_EL[res - (res >= 4095) * 4095];
        }

        return result;
    }
    
    inline FieldElement operator~ () const
    {
        FieldElement result(INV_EL[m_iElement-1]);
        return result;
    }
private:
    fe_type m_iElement;
};

const FieldElement ZERO_ELEMENT = FieldElement(0);
const FieldElement ONE_ELEMENT = FieldElement(1);

#endif FIELDELEMENT_H