float dR(float eta1, float phi1, float eta2, float phi2);
float dR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt( ( eta1 - eta2 )*( eta1 - eta2 ) + std::pow(TMath::ATan2(TMath::Sin( phi1 - phi2), TMath::Cos(phi1-phi2)),2) );
}
