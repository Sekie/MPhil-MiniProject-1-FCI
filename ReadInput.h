class InputObj
{
    public:
        void GetInputName();
        void Set();
        int aElectrons;
        int bElectrons;
        int aOrbitals;
        int bOrbitals;
        std::map< std::string, double > Integrals;
        std::string OutputName;
    private:
        std::string InputName;
};