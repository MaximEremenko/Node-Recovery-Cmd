#ifndef FIELDELEMENT_H
#define FIELDELEMENT_H

#include "FieldDefs.h"

class LogFieldElement;

class FieldElement {
public:
    FieldElement()
    {
        m_iElement = 0;
    }
        
    explicit FieldElement(fe_type iElement)
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
    
    inline const FieldElement operator= (const FieldElement element)
    {
        m_iElement = element.m_iElement;
        return *this;
    }
    
    inline const FieldElement operator= (const fe_type iElement)
    {
        m_iElement = iElement;
        return *this;
    }
    
    inline bool operator == (const FieldElement element) const
    {
        return (m_iElement == element.m_iElement);
    }
    
    inline FieldElement operator+ (const FieldElement element) const
    {
        FieldElement result(m_iElement^element.m_iElement);
        return result;
    }
    
    inline FieldElement operator+ (const fe_type iElement) const
    {
        FieldElement result(m_iElement^iElement);
        return result;
    }
    
    inline const FieldElement operator+= (const FieldElement element)
    {
        m_iElement ^= element.m_iElement;
        return *this;
    }
    
    inline const FieldElement operator+= (const fe_type iElement)
    {
        m_iElement ^= iElement;
        return *this;
    }

/*    inline FieldElement operator* (const FieldElement element) const
    {
        /*printf("fe*: x = %x, y= %x, logX = %x, logY = %x, logx + logy = %x, alog[logx + logy] = %x\n",
            m_iElement,
            element.m_iElement,
            LOG16[m_iElement],
            LOG16[element.m_iElement],
            LOG16[m_iElement] + LOG16[element.m_iElement],
            ALOG16[LOG16[m_iElement] + LOG16[element.m_iElement]]); * /
            

        FieldElement result(GF_W16_INLINE_MULT(LOG16, ALOG16, m_iElement, element.m_iElement));

        return result;
    }
    
    inline FieldElement operator* (const fe_type iElement) const
    {
        FieldElement result(GF_W16_INLINE_MULT(LOG16, ALOG16, m_iElement, iElement));

        return result;
    }
    */
    
    inline LogFieldElement toLog() const;

private:
    fe_type m_iElement;
    friend class NonZeroFieldElement;
};

const unsigned int ZeroFlag = 0xFF000000;
const unsigned int ZeroInit = 0x01000000;

unsigned int maxLog = 0;

class LogFieldElement
{

    LogFieldElement(unsigned int val) : m_iLogElement(val) {}
public:
    LogFieldElement() = default;

    explicit LogFieldElement(FieldElement elem) : m_iLogElement(LOG16_UI[elem.getElement()]) {}

    FieldElement toNormal() const
    {
        if (m_iLogElement & ZeroFlag)
        {
            return FieldElement(0);
        }
        /*printf("toNormal: logElem = %x, index = %x, normal = %x\n",
            m_iLogElement,
            m_iLogElement & 0x1FFFF,
            ALOG16[m_iLogElement & 0x1FFFF]); */
//        if ((m_iLogElement & ~ZeroFlag) > maxLog)
        //{
          //  maxLog = (m_iLogElement & ~ZeroFlag);
        //}

        return FieldElement(ALOG16[correctLog(m_iLogElement & ~ZeroFlag)]);
    }

    LogFieldElement operator * (const LogFieldElement other) const
    {
        return  m_iLogElement + other.m_iLogElement;
    }

    LogFieldElement operator * (const FieldElement other) const
    {
        return LogFieldElement(m_iLogElement + LOG16_UI[other.getElement()]);
    }

    LogFieldElement operator * (const fe_type other) const
    {
        return LogFieldElement(m_iLogElement + LOG16_UI[other]);
    }

    LogFieldElement operator / (const LogFieldElement other) const
    {
        // Предполагаем, что other != 0, так как тут так никогда не бывает.
        int result = static_cast<int>(m_iLogElement & ~ZeroFlag) - other.m_iLogElement;
        while(result < 0)
            result += 0xFFFF;
        return (m_iLogElement & ZeroFlag) | static_cast<unsigned int>(result);
    }

    LogFieldElement operator ~() const
    {
        int result = 0 - static_cast<int>(m_iLogElement);
        while (result < 0)
            result += 0xFFFF;
        return static_cast<unsigned int>(result);
    }

    static inline int correctLog(int x)
    {
//        while (x >= 0xFFFF) x -= 0xFFFF;
        return x % 0xFFFF;
        //return x >= 0xFFFF ? x - 0xFFFF : x;
        //return (x & 0xFFFF) + (x >> 16);
    }

public:
    static LogFieldElement logOneElement;


private:

    unsigned int m_iLogElement = ZeroInit;
};

LogFieldElement FieldElement::toLog() const
{
    //printf("toLog: x = %x, logx = %x\n", m_iElement, LOG16[m_iElement]);
    return LogFieldElement(*this);
}

static const FieldElement ZERO_ELEMENT = FieldElement(0);
static const FieldElement ONE_ELEMENT = FieldElement(1);

#endif FIELDELEMENT_H