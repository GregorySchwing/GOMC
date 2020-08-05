// Author Bruce Dawson
// https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/

union Float_t
{
    Float_t(float num = 0.0f) : f(num) {}
    // Portable extraction of components.
    bool Negative() const { return i < 0; }
    int32_t RawMantissa() const { return i & ((1 << 23) - 1); }
    int32_t RawExponent() const { return (i >> 23) & 0xFF; }
 
    int32_t i;
    float f;
#ifdef _DEBUG
    struct
    {   // Bitfields for exploration. Do not use in production code.
        uint32_t mantissa : 23;
        uint32_t exponent : 8;
        uint32_t sign : 1;
    } parts;
#endif
};

union Double_t
{
    Double_t(double num = 0.0) : d(num) {}
    // Portable extraction of components.
    bool Negative() const { return i < 0; }
    int64_t RawMantissa() const { return i & (((int64_t)1 << 52) - 1); }
    int64_t RawExponent() const { return (i >> 52) & 0x7FF; }
    // 0x7FF -> 11 1's
 
    int64_t i;
    double d;
#ifdef _DEBUG
    struct
    {   // Bitfields for exploration. Do not use in production code.
        uint64_t mantissa : 52;
        uint64_t exponent : 11;
        uint64_t sign : 1;
    } parts;
#endif
};
 
bool AlmostEqualUlps(float A, float B, int maxUlpsDiff)
{
    Float_t uA(A);
    Float_t uB(B);
 
    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }
 
    // Find the difference in ULPs.
    int ulpsDiff = abs(uA.i - uB.i);
    if (ulpsDiff <= maxUlpsDiff)
        return true;
 
    return false;
}

bool AlmostEqualUlps(double A, double B, int maxUlpsDiff)
{
    Double_t uA(A);
    Double_t uB(B);
 
    // Different signs means they do not match.
    if (uA.Negative() != uB.Negative())
    {
        // Check for equality to make sure +0==-0
        if (A == B)
            return true;
        return false;
    }
 
    // Find the difference in ULPs.
    int ulpsDiff = abs(uA.i - uB.i);
    if (ulpsDiff <= maxUlpsDiff)
        return true;
 
    return false;
}