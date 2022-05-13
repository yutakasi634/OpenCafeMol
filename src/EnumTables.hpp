#ifndef OPEN_AICG2_PLUS_ENUM_TABLES_HPP
#define OPEN_AICG2_PLUS_ENUM_TABLES_HPP

// define amino acid type
enum class AAType : std::uint8_t
{
    ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE,
    LEU, LYS, MET, PHE, PRO, SER, THR, TRP, TYR, VAL
};

// define functions for range based for
AAType begin(AAType) {return AAType::ALA;}
AAType end  (AAType) {return AAType::VAL;}
AAType operator* (AAType aa_type) {return aa_type;}
AAType operator++(AAType aa_type) {
    return aa_type = AAType(static_cast<std::uint8_t>(aa_type)+1);
}

#endif // OPEN_AICG2_PLUS_ENUM_TABLES_HPP
