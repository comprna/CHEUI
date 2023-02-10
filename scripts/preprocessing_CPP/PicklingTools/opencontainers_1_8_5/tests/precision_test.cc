
// Test to make sure we don't lose any precision if we prettyPrint/
// WriteValToValue for doubles.

#include "ocval.h"
#include "stdlib.h"

#if defined(OC_FORCE_NAMESPACE) 
using namespace OC;
#endif

int main ()
{
  {
    real_8 a = 1.2345678912345654321;
    Val v = a;
    
    char filename[] = "/tmp/junkreal";
    // mktemp(filename);
    WriteValToFile(v, filename);
    Val vresult;
    ReadValFromFile(filename, vresult);
    
    
    real_8 result = vresult; 
    if (a==result) {
      cout << "OK" << endl;
    } else {
      // This happens when DBL precision is 16, not 17!
      cout << "ERROR! Precision lost going to string!" << string(Val(v)) << " " << string(Val(result)) << endl;
      int_u8* vp = (int_u8*)&a;
      int_u8* rp = (int_u8*)&result;
      cout << " bits:" << *vp << " " << *rp << endl;
    }
  }
  

  // These fail for 8 places .. i.e., 8 places isn't enough!
  //values:1.0000008e-36 1.0000008e-36
  // bits: 61482029 61482030
  {
    int_u4 bbb = 61482029;
    real_4 *ap = (real_4*)&bbb;  // 1.0000008e-36;
      // 6.103515652e-005;  // fails for 7, works for 8
    //1.2345678912345654321; // fails for 7, works for 8
    Val v = *ap;
    
    char filename[] = "/tmp/junkreal2";
    // mktemp(filename);
    WriteValToFile(v, filename);
    Val vresult;
    ReadValFromFile(filename, vresult);
    
    
    real_4 result = vresult; 
    if (*ap==result) {
      cout << "OK" << endl;
  } else {
      // This happens when DBL precision is 16, not 17!
      cout << "ERROR! Precision lost going to string!" << string(Val(v)) << " " << string(Val(result)) << endl;
      int_u4* rp = (int_u4*)&result;
      cout << " bits:" << *ap << " " << *rp << endl;
    }
    
  }

  /* This is a nice search to find things that won't work with only 8 places
     for float 
  for (int_u4 ii=0; ii!=0xFFFFFFFF; ++ii) {
    real_4* rp = (real_4*)&ii;
    

    string s = string (Val(*rp));
    real_4 out = Val(s);

    if (ii % 1000000 ==0) cout << "(" << ii << "," << *rp << "," << s << ")\n";
    
    if (out != *rp) {
      cout << endl;
      cout << "values:" << string(Val(*rp)) << " " << string(Val(out)) << endl;
      cout << "bits: " << ii << " " <<  * ((int_4*)&out) << endl;
      // exit(1);
    }
  }
  */
}
