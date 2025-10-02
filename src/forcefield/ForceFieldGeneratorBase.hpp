#ifndef OPEN_AICG2_PLUS_FORCE_FIELD_GENERATOR_BASE_HPP
#define OPEN_AICG2_PLUS_FORCE_FIELD_GENERATOR_BASE_HPP

#include <OpenMM.h>

#include <memory>
#include <string>

class ForceFieldGeneratorBase
{
  public:
    virtual ~ForceFieldGeneratorBase() = default;

    virtual std::unique_ptr<OpenMM::Force> generate() const = 0;
    virtual std::string                    name()     const = 0;
};

#endif // OPEN_AICG2_PLUS_FORCE_FIELD_GENERATOR_BASE_HPP
