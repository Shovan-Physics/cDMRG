//
// Created by Erich Mueller in 2019
// Last edited by Shovan Dutta on 7/28/21
//

#ifndef _ITENSOR_cSite_
#define _ITENSOR_cSite_


namespace itensor {

class cSite;


    // New structure for array elements
    using inttriple = std::tuple <int, int, double>;
    using quadruple = std::tuple <int, int, int, double>;
    
    
class cSite
    {
    Index s;  //Private
    Index t;  //extra index for polynomial coefficients of psi(x) and density(x)
    std::vector <inttriple> Lr; //psi at left edge
    std::vector <inttriple> Rr; //psi at right edge
    std::vector <inttriple> Hr; //local Hamiltonians
    std::vector <quadruple> Pr; //psi(x)
    std::vector <quadruple> Dr; //density(x)
    int len;
    public:
        

        // Constructor
        cSite(std::vector<int> qnums,std::vector <inttriple> psiLrules,std::vector <inttriple> psiRrules,std::vector <inttriple> Hrules,std::vector<int> auxqnums,std::vector <quadruple> psirules,std::vector <quadruple> nrules,Args const& args = Args::global())
        {
            Lr=psiLrules;
            Rr=psiRrules;
            Hr=Hrules;
            Pr=psirules;
            Dr=nrules;
            auto n = 1;
            auto tags = TagSet("Site,Num");
            if(args.defined("SiteNumber") )
            {
                n = args.getInt("SiteNumber");
                tags.addTags("n="+str(n));
            }
            len=qnums.size(); //number of basis states
            auto qn = makeqn(qnums); //store quantum numbers
            s=Index(std::move(qn),tags); //site (segment) index
            auto qnaux = makeqn(auxqnums); //null qns for polynomial coefficients
            tags.replaceTags("Site","Vector"); //represents different powers of x
            t=Index(std::move(qnaux),tags); //monomial index
        }
        
        //
        // Construct operators from sparse-array rules
        //
        ITensor OpFromRules(std::vector <inttriple> rules) const
        {
            auto sP = prime(s);
            auto op=ITensor(dag(s),sP);
            auto numrules=rules.size(); //number of nonzero elements
            for(int i=0;i<=numrules-1;i++)
            {
                auto [x,y,val] = rules[i]; //assign array elements
                op.set(s(x),sP(y),val);
            }
            return op;
        }
        
        ITensor vecOpFromRules(std::vector <quadruple> rules) const
        {
            auto sP = prime(s);
            auto op=ITensor(dag(s),sP,t);
            auto numrules=rules.size();
            for(int i=0;i<=numrules-1;i++)
            {
                auto [x,y,z,val]=rules[i]; //z gives monomial degrees
                op.set(s(x),sP(y),t(z),val);
            }
            return op;
        }
        
        //
        // Assign particle numbers to each site index
        //
        Index::qnstorage makeqn(std::vector<int> qnums)
        {
            Index::qnstorage result;
            std::vector<QN> qnlist;
            std::vector<int> counts;
            len=qnums.size(); //number of basis states
            //println("Generating Quantum numbers");
            auto qn=qnums[0];
            //println("first qn");
            //println(qn);
            qnlist.push_back(QN({"N",qn})); //particle-number labels
            counts.push_back(1); //how many states with the same qn
            int ind=0; //index for distinct qns
            for(int i=1;i<=len-1;i++) //loop over basis states
            {
                auto newqn=qnums[i];
                if (newqn==qn)
                {
                    //println("repeat");
                    //println(qn);
                    counts[ind]+=1; //if same qn, increase count by 1
                }
                else
                {
                    qn=newqn;
                    //println("new");
                    //println(qn);
                    qnlist.push_back(QN({"N",qn})); //generate new label
                    counts.push_back(1); //different qn
                    ind +=1; //different index
                }
                
            }
            
            auto unique=counts.size(); //number of distinct qns
            //println("counts");
            for(int i=0;i<=unique-1;i++)
            {
                //print(i);
                //print(" ");
                //println(counts[i]);
                result.push_back(QNInt(qnlist[i],counts[i])); //distinct labels
            }
            
            return result;
        }

        
    // Default constructor -- not very useful, except for saving a memory address
    cSite()
        { }
        
        


    // Method that gives access to index
    Index
    index() const { return s; }
        
        
        // Call with state name, and it returns the IQIndex pointer
        // We will just use integers as names
        // ITensor convention has the name being a string, so we just convert
        IndexVal
        state(std::string const& state)
        {
            return s(std::stoi(state));
        }



    // Create various local operators
	ITensor
	op(std::string const& opname, Args const& args = Args::global()) const
        {
        auto sP = prime(s);

        ITensor Op(dag(s),sP);

        if(opname == "Iden")
            {
                for (auto j=1;j<=len;j+=1)
                {Op.set(s(j),sP(j),1);}
            }
        else
		if(opname == "psiL")
            {
                Op=OpFromRules(Lr);
            }
        else
        if(opname == "psiR")
            {
                Op=OpFromRules(Rr);
            }
        else
        if(opname == "psiLdag")
            {
                Op=OpFromRules(Lr);
                Op=dag(Op).swapTags("1","0");
            }
        else
        if(opname == "psiRdag")
            {
                Op=OpFromRules(Rr);
                Op=dag(Op).swapTags("1","0");
            }
        else
        if(opname == "H")
            {
                Op=OpFromRules(Hr);
            }
        else
        if(opname == "psiLdagpsiL")
            {
                auto Op1=OpFromRules(Lr);
                auto Op2=swapTags(dag(Op1),"1","0");
                Op=prime(Op2)*Op1;
                Op.swapTags("2","1");
            }
        else
        if(opname == "psiRdagpsiR")
            {
                auto Op1=OpFromRules(Rr);
                auto Op2=swapTags(dag(Op1),"1","0");
                Op=prime(Op2)*Op1;
                Op.swapTags("2","1");
            }
        else
        if(opname == "psi")
            {
                Op=vecOpFromRules(Pr);
            }
        else
        if(opname == "psidag")
            {
                Op=vecOpFromRules(Pr);
                Op=dag(Op).swapTags("Site,1","Site,0");
            }
        else
        if(opname == "n")
            {
                Op=vecOpFromRules(Dr);
            }
        else
            {
            Error("Operator " + opname + " name not recognized");
            }
        return Op;
        }
     
        
    }; //cSite class

    
    // Create site set for constructing MPS and MPOs
    // template<typename SiteType>
    class BSiteSet : public SiteSet
    {
    public:

        BSiteSet() { }


        BSiteSet(SiteStore && sites)
        {
            SiteSet::init(std::move(sites));
        }

    };
    

} //namespace itensor

#endif
