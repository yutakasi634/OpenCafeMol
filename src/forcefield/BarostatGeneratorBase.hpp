#ifndef OPEN_AICG2_PLUS_BAROSTAT_GENERATOR_BASE_HPP
#define OPEN_AICG2_PLUS_BAROSTAT_GENERATOR_BASE_HPP

#include <memory>
#include <OpenMM.h>

class BarostatGeneratorBase
{
  public:
    virtual ~BarostatGeneratorBase() = default;

    virtual std::unique_ptr<OpenMM::Force> generate()    const = 0;
    virtual double                         temperature() const = 0;
    virtual std::size_t                    frequency()   const = 0;
    virtual std::string                    name()        const = 0;
    virtual void                           dump_info()   const = 0;
};

#endif // OPEN_AICG2_PLUS_BAROSTAT_GENERATOR_BASE_HPP
