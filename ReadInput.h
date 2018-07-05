class InputObj
{
    public:
        void GetInputName();
        void SetNames(char*, char*);
        void Set();
        int aElectrons;
        int bElectrons;
        int aOrbitals;
        int bOrbitals;
        int NumberOfEV;
        std::map< std::string, double > Integrals;
        std::string OutputName;
        std::string InputName;
		std::string TruncatedCI = "FCI";
};
