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
    
    fe_type getElement() const
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
        //FieldElement result;

        //FieldElement result(GF.multiply.w32(&GF, m_iElement, element.m_iElement));

        /*printf("fe*: x = %x, y= %x, logX = %x, logY = %x, logx + logy = %x, alog[logx + logy] = %x\n",
            m_iElement,
            element.m_iElement,
            LOG16[m_iElement],
            LOG16[element.m_iElement],
            LOG16[m_iElement] + LOG16[element.m_iElement],
            ALOG16[LOG16[m_iElement] + LOG16[element.m_iElement]]);*/
            

        FieldElement result(GF_W16_INLINE_MULT(LOG16, ALOG16, m_iElement, element.m_iElement));

        //if (m_iElement && element.m_iElement)
        //{
            //fe_type a_deg, b_deg, res = 0;
            //a_deg = EL_TO_DEG[m_iElement - 1];
            //b_deg = EL_TO_DEG[element.m_iElement - 1];
            //res = (a_deg + b_deg);
            //result.m_iElement = DEG_TO_EL[res - (res >= 4095)*4095];
            //result.m_iElement = DEG_TO_EL[EL_TO_DEG[m_iElement - 1] + EL_TO_DEG[element.m_iElement - 1]];
        //}

        return result;
    }
    
    inline FieldElement operator* (const fe_type iElement) const
    {
        //FieldElement result;

        //FieldElement result(GF.multiply.w32(&GF, m_iElement, iElement));
        FieldElement result(GF_W16_INLINE_MULT(LOG16, ALOG16, m_iElement, iElement));

        //if (m_iElement && iElement)
        //{
            //fe_type a_deg, b_deg, res = 0;
            //a_deg = EL_TO_DEG[m_iElement - 1];
            //b_deg = EL_TO_DEG[iElement - 1];
            //res = (a_deg + b_deg);
            //result.m_iElement = DEG_TO_EL[res - (res >= 4095) * 4095];
            //result.m_iElement = DEG_TO_EL[EL_TO_DEG[m_iElement - 1] + EL_TO_DEG[iElement - 1]];
        //}

        return result;
    }
    
    inline FieldElement operator~ () const
    {
        //FieldElement result(INV_EL[m_iElement-1]);
        //FieldElement result(GF.inverse.w32(&GF, m_iElement));
        FieldElement result(GF_W16_INLINE_DIV(LOG16, DALOG16, 1, m_iElement));
        return result;
    }
private:
    fe_type m_iElement;
    friend class NonZeroFieldElement;
};

#include "iostream"

void assert_x(fe_type iElement)
{
    if (1 / iElement)
        std::cout << "";
}

class LogFieldElement;

class NonZeroFieldElement {
public:
    NonZeroFieldElement(fe_type iElement)
    {
        m_iElement = iElement;
        //assert_x(m_iElement);
    }

    NonZeroFieldElement(const FieldElement& element)
    {
        m_iElement = element.m_iElement;
        //assert_x(m_iElement);
    }

    NonZeroFieldElement(const  NonZeroFieldElement& element) = default;

    ~NonZeroFieldElement()
    {
    }

    FieldElement toFieldElement()
    {
        return FieldElement(getElement());
    }

    fe_type getElement() const
    {
        return m_iElement;
    }

    inline const NonZeroFieldElement& operator= (const FieldElement& element)
    {
        m_iElement = element.m_iElement;
        //assert_x(m_iElement);
        return *this;
    }

    inline const NonZeroFieldElement& operator= (const fe_type iElement)
    {
        m_iElement = iElement;
        //assert_x(m_iElement);
        return *this;
    }

    inline bool operator == (const NonZeroFieldElement& element) const
    {
        return (m_iElement == element.m_iElement);
    }

    inline FieldElement operator+ (const NonZeroFieldElement& element) const
    {
        FieldElement result(m_iElement ^ element.m_iElement);

        return result;
    }

    inline FieldElement operator+ (const fe_type iElement) const
    {
        FieldElement result(m_iElement ^ iElement);
        return result;
    }


    inline LogFieldElement toLog() const;
    

    inline NonZeroFieldElement operator~ () const
    {
        //FieldElement result(INV_EL[m_iElement-1]);
        //FieldElement result(GF.inverse.w32(&GF, m_iElement));
        NonZeroFieldElement result(GF_W16_INLINE_DIV(LOG16, DALOG16, 1, m_iElement));
        return result;
    }
private:
    fe_type m_iElement;
};


class LogFieldElement
{
public:
    LogFieldElement(int iLogElement) : m_iLogElement(iLogElement) {}

    NonZeroFieldElement toNormal()
    {
/*        printf("toNormal: logElem = %x, index = %x, normal = %x\n",
            m_iLogElement,
            m_iLogElement & 0x1FFFF,
            ALOG16[m_iLogElement & 0x1FFFF]); */
        return NonZeroFieldElement(ALOG16[correctLog(m_iLogElement)]);
    }

    LogFieldElement operator * (const LogFieldElement& other)
    {
        /*printf("log*: logx = %x, logy = %x, result = %x\n",
            m_iLogElement,
            other.m_iLogElement, 
            LogFieldElement(m_iLogElement + other.m_iLogElement).m_iLogElement);*/
           
        return LogFieldElement((m_iLogElement + other.m_iLogElement));
    }

    FieldElement operator * (const FieldElement& other)
    {
        return FieldElement(ALOG16[(m_iLogElement + LOG16[other.getElement()])]);
    }

    NonZeroFieldElement operator * (const NonZeroFieldElement& other)
    {
        return NonZeroFieldElement(ALOG16[(m_iLogElement + LOG16[other.getElement()])]);
    }

    static inline int correctLog(int x)
    {
        //return x >= 0xFFFF ? x - 0xFFFF : x;
        return (x & 0xFFFF) + (x >> 16);
    }

private:

    int m_iLogElement;
};

LogFieldElement NonZeroFieldElement::toLog() const
{
    //printf("toLog: x = %x, logx = %x\n", m_iElement, LOG16[m_iElement]);
    return LogFieldElement(LOG16[m_iElement]);
}

FieldElement operator + (const FieldElement& left, const NonZeroFieldElement& right)
{
    return FieldElement(left.getElement() ^ right.getElement());
}

FieldElement operator + (const NonZeroFieldElement& left, const FieldElement& right)
{
    return FieldElement(left.getElement() ^ right.getElement());
}

static const FieldElement ZERO_ELEMENT = FieldElement(0);
static const FieldElement ONE_ELEMENT = FieldElement(1);

#endif FIELDELEMENT_H