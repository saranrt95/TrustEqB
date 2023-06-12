# include <iostream>
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <list>
# include <iomanip>
#include <queue>
#include <time.h>
#include "rulesel.c"

////////////////////////////// Grandezze di #define base /////////////////////////////////////////////////
#define UnKilobit 1000
#define UnMegabit UnKilobit * 1000
#define UnByte 8 //in bit (ovviamente)
#define DimElementare UnByte //in bit (ovviamente)
#define CellOverHead 4 // in byte
#define CellSize 188 // in byte

const double DimPacchetto= (80*UnByte);  //in bit per il GeneraConBurst

#define MaxNumConn_IP 121

int MaxNumConn_IP_AttiveOra;

const int NumBuffers_IP=2;

//Vincoli di SLA fissati poi in Inizializza_IP
double ReferenceIdealeLoss_0; 
double ReferenceIdealeLoss_1; 

double MaxDelay[NumBuffers_IP];

double ReferenceIdealeDelay_0; 
double ReferenceIdealeDelay_1; 
//END Vincoli di SLA fissati poi in Inizializza_IP

///////////////////////////////////////////////////////////////////// END Grandezze di #define base /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// Per la statistica pacchetto IP /////////////////////////////////////////////////////////////////////

//per generare pack con dim=1byte mettere:  DimensioneMediaPacchetto=1, MaxDimensionePacchetto=1*DimensioneMediaPacchetto
#define DimensioneMediaPacchetto 3 // IN BYTE!!!
#define MaxDimensionePacchetto 2*DimensioneMediaPacchetto
#define DimDistrBimodaleMIN 40 
#define DimDistrBimodaleMAX 540

// Distribuzione Trimodale (a,b,c,pa,pb); papero Ajmone10/02 Trans. on Net. TRIMODAL(47 byte, 576 byte, 1500 byte; 0.559, 0.200)
#define a_DimDistrTrimodale 47
#define b_DimDistrTrimodale 576
#define c_DimDistrTrimodale 1500
double pa_DistrTrimodale = 0.559;
double pb_DistrTrimodale = 0.200;
double a_PacchettiArrivatiDistrTrimodale[MaxNumConn_IP];
double b_PacchettiArrivatiDistrTrimodale[MaxNumConn_IP];
double c_PacchettiArrivatiDistrTrimodale[MaxNumConn_IP];
//END Distribuzione Trimodale

//Per controllare la media della distribuzione uniforme
int ContDimPackArrDistrUniforme[MaxNumConn_IP][MaxDimensionePacchetto+1];
double DimMediaEffettivaDistrUniforme[MaxNumConn_IP];
//END per controllare la media della distribuzione uniforme

//Per controllare la media della distribuzione esponenziale
double DimMediaEffettivaDistrEXP[MaxNumConn_IP];
//END per controllare la media della distribuzione esponenziale

//Per controllare la media della distribuzione Bimodale
double DimMediaEffettivaDistrBimodale[MaxNumConn_IP];
int PacchettiPiccoliBimodale[MaxNumConn_IP], PacchettiGrossiBimodale[MaxNumConn_IP];
double SogliaXDistrBimodale=0.08;
//END per controllare la media della distribuzione Bimodale

///////////////////////////////////////////////////////////////////// END Per la statistica pacchetto IP /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// Lato IP /////////////////////////////////////////////////////////////////////

int MaxDimBuffer_IP[NumBuffers_IP];      //IN BYTE!!
int DimBufferAttuale_IP[NumBuffers_IP];//IN BYTE!!
bool Buffer_IPVuoto[NumBuffers_IP];

/*attenzione sarebbe molto interessante avere una analisi dei successivi 3 contatori anche per dimpacchetto ossia
PacchettiX[MaxNumConn_IP][MaxDimensionePacchetto] X=Arrivati, Persi, Serviti*/
double PacchettiArrivati[MaxNumConn_IP];
double PacchettiPersi[MaxNumConn_IP];

double PacchettiArrivatiTot[NumBuffers_IP];
double PacchettiPersiTot[NumBuffers_IP];
double PacchettiServitiTot[NumBuffers_IP];
double Bit_PacchettiServitiTot[NumBuffers_IP];
double PacchettiServitioverATMTot[NumBuffers_IP];

//double Ritardo_IP[MaxNumConn_IP];
//double Ritardo_IPTot;

double PPerdita_IP[NumBuffers_IP];
double DelayPackAttuale;
double UltimoDelayFlussoIPoATM[NumBuffers_IP];

double PPerdita_IP_IDEALE; // da passare alla riallocazione con banda equivalente bufferless

//double PPerdita_IPSingolaConn[MaxNumConn_IP];
//double DiffRispettoPPerdita_IP[MaxNumConn_IP];

double PacchettiTotaliPersiNellaTramaATM;
double PacchettiTotaliPersiNellaTramaATM_FlussoIP[NumBuffers_IP];
double PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[NumBuffers_IP];

// per fissare la statistica sulla dimensione del pacchetto
bool DistrUniformePack=false;
bool DistrExpPack=false;
bool DistrBimodalePack=false;
bool DistrTrimodalePack=true;
//END per fissare la statistica sulla dimensione del pacchetto

//riguardo a validazione misura PPerdita_IP
bool StatisticaPPerdita_IPvalidata=false;
double MaxScostamentoPPerd;
double DevStandPopPPerd;
//END riguardo a validazione misura PPerdita_IP

// Per GeneraConBurst il traffico IP
bool GeneraConBurst=true;
double MediaDurataSilenzio[NumBuffers_IP]; 
double MediaDurataBurst[NumBuffers_IP]; 
double Bp[MaxNumConn_IP];
double Init_Bp[NumBuffers_IP];

//CLAUDIO --> pareto
bool GeneraSecondoClaudio=0;
double alfa;
double deltaMediaDurataBurst[NumBuffers_IP];
double deltaMediaDurataSilenzio[NumBuffers_IP];
//END CLAUDIO
//MARIO --> MMDP
bool GeneraSecondoMario=1;
double MediaDurataBurstInCelle;
//END MARIO

double DurataAttualeBurst[MaxNumConn_IP];
double IstanteFineAttualeBurst[MaxNumConn_IP];

int NumConnAttive;
double NumBurstGenerati;
double DurataEffettivaDeiBurst;
double NumSilenziGenerati;
double DurataEffettivaDeiSilenzi;
// END Per GeneraConBurst il traffico IP

//Generazione video

FILE *FileTracciatoVideo_Pacchetti, *FileTracciatoVideo_Tempi;

const int NumCampioniTracciatoVideo=	25407	;
double CampioniTracciatoVideo_Pacchetti[NumCampioniTracciatoVideo];
double CampioniTracciatoVideo_Tempi[NumCampioniTracciatoVideo];
int IndiceTracciatoVideo=0;

//END Generazione video

///////////////////////////////////////////////////////////////////// END Lato IP /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// Lato ATM /////////////////////////////////////////////////////////////////////

int DimBufferAttuale_ATM;
int DimBufferAttuale_ATM_scelto;
double CelleArrivate;
double CellePerse;
double CelleServite;

double BitArrivati_monitoraggioATM;
double CampioneBitRate_monitoraggioATM;
double CampioneBitRateQuadro_monitoraggioATM;

int MaxDimBuffer_ATM;

double ScartoQuadraticoMedio_IPTramaATM=0.0;
//devo conoscere le dimensioni medie x i conti di RCBC successivi
double DimPacchettoMedio[NumBuffers_IP];

double LossRate_IDEALI_0;
double LossRate_IDEALI_1;
double DelayRate_IDEALI_0;
double DelayRate_IDEALI_1;

double UltimoIstantePerdita_ATM;
double UltimoIstantePerdita_IPoATM[NumBuffers_IP];
double UltimoIstantePerdita_IPoATM_Delay[NumBuffers_IP];

bool Buffer_ATMVuoto;

bool BusyPeriodAttivo[NumBuffers_IP];
double InizioBusyPeriod_ATM;
double InizioBusyPeriod_IPoATM[NumBuffers_IP];

double ContributoIPA_ATM;
double ContributoIPA_IPoATM[NumBuffers_IP];
double ContributoIPA_IPoATM_Delay[NumBuffers_IP];

double StimaIPADerivata_ATM;
double StimaIPADerivata_IPoATM[NumBuffers_IP];
double StimaIPADerivata_IPoATM_Delay[NumBuffers_IP];

bool PackInCorsoTramaATMGiaPerso[MaxNumConn_IP];
double PacchettiPersiNellaTramaATM[MaxNumConn_IP];
double PacchettiInRitardoNellaTramaATM[MaxNumConn_IP];

double PacketLossTramaATM;
double PacketLossTramaATM_FlussoIP[NumBuffers_IP];
double TotSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[NumBuffers_IP];
double TotQuadroSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[NumBuffers_IP];
double ScostamentoPercentualeLoss_IPoATM[NumBuffers_IP];
double TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[NumBuffers_IP];
double PacketDelayTramaATM_FlussoIP[NumBuffers_IP];
double ScostamentoPercentualePdelay_IPoATM[NumBuffers_IP];

bool AggiornamentoPerArrivoInConvergenza=false;
bool NonEntrarePiu=false;
double TempoPerArrivareInConvergenza=-1.0;//metto -1.0 appositamente per accorgermi se succede che non riesce a calcolarlo
double MediaPrimaConvergenzaScostamentoPercentualeLoss_IPoATM[NumBuffers_IP];
int EntrateInRiallocazione=0;
int EntrateInRiallocazioneIdealDelay=0;

int NumVolteSopraTarget[NumBuffers_IP];
double MediaDifferenzaSopraTarget[NumBuffers_IP];

///////////////////////////////////////////////////////////////////// END Lato ATM /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// Variabili di controllo Relay Entity /////////////////////////////////////////////////////////////////////

//Ocio che Ct_IP puo' essere inizializzato anche in Inizializza_IP nel caso in cui avvegna un cambio di statistiche (e.g., #connessioni) durante la simulazione
double Ct_IP[NumBuffers_IP];

double CtInSorgente_IP=1.0*UnMegabit;
//double T_out_IP = ( (double)(DimElementare) ) / Ct_IP[0]; //tempo per emettere sul canale IP una cella elementare
//double T_in_IP = ( (double)(DimElementare) ) / (CtInSorgente_IP); //tempo che la sorgente impiega per inserire sul canale una cella elementare
double T_out_IP[NumBuffers_IP];
double T_in_IP[MaxNumConn_IP];

double Ct_ATM=0.0;
double Ct_ATM__OLD=0.0;
double Ct_ATM_Tot=0.0;
double Ct_ATM_TotQuadro=0.0;
double Ct_ATM_IPA[NumBuffers_IP];
double T_out_ATM =  (((double)(CellSize))*(double)(DimElementare)) / Ct_ATM;   //tempo per emettere sul canale ATM una cella ATM

bool BuffersInparallelo = 0;
double T_in_ATM[NumBuffers_IP];

double PPerdita_ATM;

double DimIntervalloRiallocazioniBanda=		1*60	;

double DurataSimulazione = 3*60 ;

double NumCampioni_bitrate=10;
double DimIntervalloMonitor_bitrate=(1*DurataSimulazione/3) / NumCampioni_bitrate;
double Addestramento_Media_bitRate_ATM;
double Addestramento_StDev_bitRate_ATM;

// 0: RCBC,			1: EqB,			2: PID,			3: Ideal			4: Bm			5: IdealeLoss_2				6: IdealeDelay
int TecnicaRiallocazione=0;
bool LavoraSuPDelay=0;

int contazeri=0;

int IndiceletturaPacchettiIPServiti=0;
const int NumCampioniPacchettiIPServiti=177;//leggendo il numero di righe dei file corrispondenti
double Bit_PacchettiIPServiti_0[NumCampioniPacchettiIPServiti];
double Bit_PacchettiIPServiti_1[NumCampioniPacchettiIPServiti];
double PlossoverATM_0[NumCampioniPacchettiIPServiti];
double PlossoverATM_1[NumCampioniPacchettiIPServiti];

double CellTaxTeorico=0.0;
double Ct_ATMSecondoCellTaxTeorico=0.0;
double ProvaIncremetoCt_ATM=1.0*UnMegabit;

FILE *Time;
FILE *PlossNellaTramaATM;
FILE *PlossNellaTramaATM_0;
FILE *PlossNellaTramaATM_1;
FILE *PPerditaFlussoIP_0;
FILE *PPerditaFlussoIP_1;
FILE *ScartoQuadraticoMedio;
FILE *ValoreSpintaGradiente;
FILE *ValoreCt_ATM;
FILE *SimulatedCellTax;
FILE *ValoreCt_IP_0;
FILE *ValoreCt_IP_1;
FILE *Bit_PackIPServiti_0;
FILE *Bit_PackIPServiti_1;

///////////////////////////////////////////////////////////////////// END Variabili di controllo Relay Entity /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// Struttura simulatore /////////////////////////////////////////////////////////////////////

using namespace std ;

enum event {InviaPacchetto_NelCanale, Arrivo_Pacchetto, Arrivo_Cella_ATM,
InviaCella_NelCanale_ATM, Riallocazione_Ct_ATM, Cambio_Statistiche, AzzeraContatori, Inizio_Burst, Monitor_bitrate};

struct track
{
  double istante;
  event evento;
  int Connessione;
  int Buffer_IP;

  //se l'evento e' Arrivo_Cella_ATM
  int ConnIPGenerante;
  bool PrimaCellaDelPacchetto;
  bool UltimaCellaDelPacchetto;

  public:
  bool operator < (track); 
};

list<track>L(0);

// per conservare lo "stato" del buffer IP
struct PacchettoInserito{
	int NumConn;      //numero della connesione che lo ha generato
	int Dimensione; //in byte
    double TInserimentoNelBuffer;
};
typedef queue <PacchettoInserito> Coda;
Coda CodaPacchettiInseriti[NumBuffers_IP];
// END per conservare lo "stato" del buffer IP

// per conservare lo "stato" del buffer ATM
struct CellaInserita{
	int NumConn;      //numero della connesione che l'ha generata
	//int Dimensione; //in byte
    double TInserimentoNelBuffer;
};
typedef queue <CellaInserita> Coda2;
Coda2 CodaCelleInserite;
// END per conservare lo "stato" del buffer ATM 

///////////////////////////////////////////////////////////////////// END Struttura simulatore /////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////// PID ///////////////////////////////////////////////////////////
//tapullo del PID, in spagnolo "parche"

//per il PID
double Kp[NumBuffers_IP];
double Kd[NumBuffers_IP];
double Ki[NumBuffers_IP];

double e_k[NumBuffers_IP];
double e_kMeno_1[NumBuffers_IP];
double e_kMeno_2[NumBuffers_IP];

double Ct_ATMSeparati[NumBuffers_IP];
/////////////////////////////////////////////////////////// END PID ///////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// FUNZIONI /////////////////////////////////////////////////////////////////////

//Lato IP
void ArrivoPacchetto(int NumConn, double ClockAttuale);
void InviaPacchettoNelCanale(double ClockAttuale, int Buffer_IP);
int CheBufferIPCorrente(int NumConn);
double DistribuzioneUniforme0_max(int max);
void Inizializza_IP();
void AggiornaContatoriTot_IP();
void CalcolaPPerdita_IP();
void StampaRisultati_IP();
int GenAttualeDimPack(int NumConn);
void InizioBurst(int NumConn, double ClockAttuale);
double GenDurataBurst(int NumConn);
double GenDurataSilenzio(int NumConn);
void InizializzaTracciatoVideo();
void InizializzaTracciatoPacchettiIPServiti();
void InizializzaTracciatoPloss();

//Lato ATM
void Inizializza_ATM();
void InizializzaT_in_ATM();
void AggiornaContatoriTot_ATM();
void GeneraCelleDaPacchetto(double ClockAttuale, int DimPacchettoAttuale, int NumConn);
void ArrivoCella_ATM(double ClockAttuale, int ConnIPGenerante, bool PrimaCellaDelPacchetto, bool UltimaCellaDelPacchetto);
void InviaCellaNelCanale_ATM(double ClockAttuale);
void CalcolaPPerdita_ATM();
void StampaRisultati_ATM();

//Utilita'
void SkipSobol(FILE** FileSobol_ref, int skip);
void Inizializza();
void StampaRisultati();
double TrovaMinimo(double *Array, int MaxDim);
int TrovaIndiceMinimo(double *Array, int MaxDim);
double TrovaMassimo(double *Array, int MaxDim);
double NumeroRand(int max);
double modulo(double num);
void AzzeraContatoriPerRiprendereIlContoDalRegime();
void Svuota_EventiRimastiInCoda();
void StampaLoss();
void RiallocazioneCt_ATM(double ClockAttuale);
void RiallocazioneCt_ATM_IPA(double ClockAttuale);
void RiallocazioneCt_ATM_BEqBufferless(double ClockAttuale);
void RiallocazioneCt_ATM_PID(double ClockAttuale);
void RiallocazioneCt_ATM_Ideale(double ClockAttuale);
void RiallocazioneCt_ATM_Bm(double ClockAttuale);
void RiallocazioneCt_ATM_IdealeLoss_2(double ClockAttuale);
void StampaStato_RiallocazioneCt_ATM(double ClockAttuale, double ScartoQuadraticoMedio_IPTramaATM, double SpintaGradiente, double Ct_ATM);
void CambioStatistiche(double ClockAttuale);
int ConvertiAdInteroSuperiore(double num);
void Monitorbitrate(double ClockAttuale);
void CalcolaMediaVarianza();
double RegoleRulex(int N, double bw, double loss, double mean, double var, int B, int MaxB, double Silence, double Bp);

///////////////////////////////////////////////////////////////////// END FUNZIONI /////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////// main() /////////////////////////////////////////////////////////////////////

void main()
{
event EventoAttuale;
int i,j,scelta,Ripetizioni, skip_sobol = 0;
double sommaprestazioni=0.0;
double Banda_daregistrare, loss_daregistrare;
double Ct_ATM_passo0;
  
srand(time(NULL));
double clock_sim;

const int TentativiDelta=3;
Ripetizioni = 1000;
const int NumInputSobol=6;
double Prestazioni[TentativiDelta];
double RandSobol[NumInputSobol];
bool supervised = false, sobol = true;
char tmp[120], output_file[120], buf[5];
const char *confidence;
const char *error = NULL;

FILE *FileSobol=fopen("sobol6x10000.txt","r");
FILE *Output;

if (supervised) {
    strcpy(output_file, "C:\\Users\\AndreaScorzoni\\Desktop\\Andrea\\Safe.Coop\\Buffer\\sw su RCBC Andrea\\Iteration\\OutTrain");
    for (i = 1; ; i++) {
      strcpy(tmp, output_file);
      sprintf(buf, "%d", i);
      strcat(tmp, buf);
      strcat(tmp, ".txt");
      if (!fopen(tmp, "r")) {
        break;
      }
    }
    strcpy(output_file, tmp);
  }
  else
    strcpy(output_file, "OutTrain.txt");

Output = fopen(output_file,"w");

fprintf(Output,"N\tmean\tvar\tloss\tbw\tB\tMaxB\tSilence\tBp\tapprox_mean\tout\n");

//per calcolo dei tempi di generazione del file di addestramento
clock_t inizio, fine;	
unsigned int durata = 0;
unsigned int durata_sec = 0;
unsigned int durata_min = 0;
unsigned int durata_ore = 0;	

SkipSobol(&FileSobol, skip_sobol);

inizio = clock();
double b, approx_mean;

for (j = 0;j < Ripetizioni;j++) {//ripetizioni

  if (sobol) {
    //estrazioni con Sobol
    for (i = 0; i < 6; i++) {
      if (i < 5)
        fscanf(FileSobol, "%lf\t", &RandSobol[i]);
      else
        fscanf(FileSobol, "%lf\n", &RandSobol[i]);
    }
    //MaxNumConn_IP_AttiveOra=69+(int)(RandSobol[0]*40)+1;	
    MaxNumConn_IP_AttiveOra = 69 + (int)(RandSobol[0] * 40) + 1;
    MaxDimBuffer_ATM = (int)(RandSobol[2] * 500) + 1;
    DimBufferAttuale_ATM_scelto = (int)(RandSobol[3] * MaxDimBuffer_ATM);
    MediaDurataSilenzio[0] = RandSobol[4] * 5.0;
    Init_Bp[0] = (10.0 + RandSobol[5] * 40.0)*UnKilobit;

    //Euristica: do la banda tra la Bm e la Bp
    b = (1 / (1.0 / (1.0 + MediaDurataSilenzio[0])));
    Ct_ATM = ((double)MaxNumConn_IP_AttiveOra) * ((Init_Bp[0] / b) + (RandSobol[1] * (Init_Bp[0] - (Init_Bp[0] / b))));
  }
  else {
    MaxNumConn_IP_AttiveOra = 69 + (int)(NumeroRand(40)) + 1;
    MaxDimBuffer_ATM = (int)(NumeroRand(500)) + 1;
    DimBufferAttuale_ATM_scelto = (int)(NumeroRand(MaxDimBuffer_ATM));
    MediaDurataSilenzio[0] = NumeroRand(5.0);
    Init_Bp[0] = (10.0 + NumeroRand(40.0))*UnKilobit;

    //Euristica: do la banda tra la Bm e la Bp
    b = (1 / (1.0 / (1.0 + MediaDurataSilenzio[0])));
    Ct_ATM = ((double)MaxNumConn_IP_AttiveOra) * ((Init_Bp[0] / b) + (NumeroRand(Init_Bp[0] - (Init_Bp[0] / b))));
  }
  approx_mean = (Init_Bp[0] / b) / (UnKilobit);

  
  if (supervised) {
    error = ApplyRules(Ct_ATM/(UnMegabit), MediaDurataSilenzio[0], Init_Bp[0] / (UnKilobit), approx_mean);
    //printf("confidence = %s", confidence);
  }
  if (supervised && !strcmp(error, "correct")) {
    j--;
    continue;
  }

  printf("\n\t#%d\n", j);
  /*for(i=0; i<6; i++)
    printf("%f",RandSobol[i]);*/

  for (i = 0; i < TentativiDelta; i++) {//ciclo su +- delta

    Inizializza();

    //Per abilitare le regole del Rulex parto da -1 nel for precedente in modo da fare un primo giro in cui calcolo media e varianza
    if (i == 0) {
      //Se uso le regole del Rulex gli passo la banda come settata nel primo (o precedente?) training	; e.g., tra min e max senza euristica, tra Bm e Bp con l'euristica BmBp					
      /*Ct_ATM=RegoleRulex(MaxNumConn_IP_AttiveOra, (Ct_ATM/(UnMegabit)),PacketLossTramaATM_FlussoIP[0],
                      (Addestramento_Media_bitRate_ATM/(UnMegabit)), (Addestramento_StDev_bitRate_ATM/(UnMegabit)),DimBufferAttuale_ATM_scelto,
                      MaxDimBuffer_ATM, MediaDurataSilenzio[0], (Init_Bp[0]/(UnKilobit)));
      Ct_ATM*=UnMegabit;*/

      /*Partendo da -1, già che calcolo la media, faccio un addestramento passandogli la media
      Ct_ATM=Addestramento_Media_bitRate_ATM;*/

      Ct_ATM_passo0 = Ct_ATM;
    }
    else if (i != -1) {
      Ct_ATM = Ct_ATM_passo0;
      if (i == 1)
        Ct_ATM += (Ct_ATM_passo0*0.15);
      if (i == 2)
        Ct_ATM -= (Ct_ATM_passo0*0.15);
    }

    T_out_ATM = ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;

    clock_sim = 0.0;
    while (clock_sim <= DurataSimulazione) {

      EventoAttuale = L.front().evento;
      clock_sim = L.front().istante;

      switch (EventoAttuale)
      {
      case Inizio_Burst: InizioBurst(L.front().Connessione, clock_sim);break;

      case InviaPacchetto_NelCanale: InviaPacchettoNelCanale(clock_sim, L.front().Buffer_IP);break;
      case Arrivo_Pacchetto: ArrivoPacchetto(L.front().Connessione, clock_sim);break;

      case Arrivo_Cella_ATM: ArrivoCella_ATM(clock_sim, L.front().ConnIPGenerante,
        L.front().PrimaCellaDelPacchetto, L.front().UltimaCellaDelPacchetto);break;
      case InviaCella_NelCanale_ATM: InviaCellaNelCanale_ATM(clock_sim);break;
      case Riallocazione_Ct_ATM: RiallocazioneCt_ATM(clock_sim);break;
      case Cambio_Statistiche: CambioStatistiche(clock_sim);break;
      case AzzeraContatori: AzzeraContatoriPerRiprendereIlContoDalRegime();break;
      case Monitor_bitrate: Monitorbitrate(clock_sim);break;
      }
      L.pop_front();

    }//END while(clock<=DurataSimulazione)

    //printf("\n\n\t***Condizione Uscita Simulazione soddisfatta ESCO***\n");
    //StampaRisultati();

    AggiornaContatoriTot_IP();
    AggiornaContatoriTot_ATM();
    CalcolaPPerdita_IP();
    CalcolaPPerdita_ATM();

    CalcolaMediaVarianza();

    printf("\nN=%d, bw=%f, loss=%f, m=%f, std=%f, B=%d, MaxB=%d, fi=%f, Bp=%f", MaxNumConn_IP_AttiveOra, (Ct_ATM / (UnMegabit)),
      PacketLossTramaATM_FlussoIP[0], (Addestramento_Media_bitRate_ATM / (UnMegabit)),
      (Addestramento_StDev_bitRate_ATM / (UnMegabit)), DimBufferAttuale_ATM_scelto,
      MaxDimBuffer_ATM, MediaDurataSilenzio[0], (Init_Bp[0] / (UnKilobit)));

    Prestazioni[i] = PacketLossTramaATM_FlussoIP[0];

    Svuota_EventiRimastiInCoda();

    if (i == 0) {
      Banda_daregistrare = Ct_ATM;
      loss_daregistrare = PacketLossTramaATM_FlussoIP[0];
    }

  }//END for for(int i=0; i<3; i++){//ciclo su +- delta

  sommaprestazioni = 0.0;
  for (i = 0; i < TentativiDelta; i++) {
    sommaprestazioni += Prestazioni[i];
    Prestazioni[i] = modulo(Prestazioni[i] - ReferenceIdealeLoss_0);
    //Prestazioni[i]=pow( (Prestazioni[i]-ReferenceIdealeLoss_0),2 );		
  }

  scelta = TrovaIndiceMinimo(Prestazioni, TentativiDelta);
  //tolleranza sull'oscillazione intorno al target: nella tolleranza scelgo di tenere la banda ferma
  /*if(Prestazioni[scelta]<(ReferenceIdealeLoss_0/2))
    scelta=0;*/
  if (sommaprestazioni == 0.0)
    scelta = 2;
  if (sommaprestazioni > (0.95 * 3))
    scelta = 1;
  printf("\n\t\tscelta=%d\n", scelta);
  
  //printf("\nCt_ATM = %f\n", Ct_ATM);

  if (scelta != 0)
    fprintf(Output, "%d\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\t%d\n", MaxNumConn_IP_AttiveOra, (Addestramento_Media_bitRate_ATM / (UnMegabit)),
    (Addestramento_StDev_bitRate_ATM / (UnMegabit)), loss_daregistrare, (Banda_daregistrare / (UnMegabit)),
      DimBufferAttuale_ATM_scelto, MaxDimBuffer_ATM, MediaDurataSilenzio[0], (Init_Bp[0] / (UnKilobit)),
      approx_mean, scelta);
  else
    j--;

}//END for(int j=0;j<10;j++){//ripetizioni

///////////////////////////////
fine = clock();	
durata = (unsigned int) (fine-inizio)/CLOCKS_PER_SEC;
durata_ore = durata / 3600;
durata_min = (durata % 3600)/60;
durata_sec = (durata % 3600)%60;
cout<<"\n\nDURATA: "<<durata_ore<<" ore\t"<<durata_min<<" minuti\t"<<durata_sec<<" secondi"<<endl<<endl;
inizio = clock();
///////////////////////////////	

//chiusura file
/*fclose(Time);
fclose(PlossNellaTramaATM);
fclose(PlossNellaTramaATM_0);
fclose(PlossNellaTramaATM_1);
fclose(PPerditaFlussoIP_0);
fclose(PPerditaFlussoIP_1);
fclose(ScartoQuadraticoMedio);
fclose(ValoreSpintaGradiente);
fclose(ValoreCt_ATM);
fclose(SimulatedCellTax);
fclose(ValoreCt_IP_0);
fclose(ValoreCt_IP_1);
if(TecnicaRiallocazione==0 || TecnicaRiallocazione==4){
	fclose(Bit_PackIPServiti_0);
	fclose(Bit_PackIPServiti_1);
}*/
//END chiusura file
fclose(FileSobol);
fclose(Output);

printf("\n");
system("pause");

} 

///////////////////////////////////////////////////////////////////// END main() /////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Inizializza(){

StatisticaPPerdita_IPvalidata=false;
MaxScostamentoPPerd=0.0;
DevStandPopPPerd=0.0;

Inizializza_IP();

Inizializza_ATM();

track TracciaAttuale;

double TempoXIlRegime=2 * DurataSimulazione / 3;

//programmo un AzzeraContatori
TracciaAttuale.evento = AzzeraContatori;
TracciaAttuale.istante = TempoXIlRegime;
L.push_back(TracciaAttuale);

//programmo il primo monitoraggio della bitrate
TracciaAttuale.evento = Monitor_bitrate;
TracciaAttuale.istante = TempoXIlRegime+DimIntervalloMonitor_bitrate;
L.push_back(TracciaAttuale);

L.sort();

	/*Time=fopen("Time.txt","w");
	PlossNellaTramaATM=fopen("PlossNellaTramaATM.txt","w");
	PlossNellaTramaATM_0=fopen("PlossNellaTramaATM_0.txt","w");
	PlossNellaTramaATM_1=fopen("PlossNellaTramaATM_1.txt","w");
	PPerditaFlussoIP_0=fopen("PPerditaFlussoIP_0.txt","w");
	PPerditaFlussoIP_1=fopen("PPerditaFlussoIP_1.txt","w");
	ScartoQuadraticoMedio=fopen("ScartoQuadraticoMedio.txt","w");
	ValoreSpintaGradiente=fopen("ValoreSpintaGradiente.txt","w");
	ValoreCt_ATM=fopen("ValoreCt_ATM.txt","w");
	SimulatedCellTax=fopen("SimulatedCellTax.txt","w");
	ValoreCt_IP_0=fopen("ValoreCt_IP_0.txt","w");
	ValoreCt_IP_1=fopen("ValoreCt_IP_1.txt","w");
	if(TecnicaRiallocazione==0 || TecnicaRiallocazione==4){
		Bit_PackIPServiti_0=fopen("Bit_PackIPServiti_0.txt","w");
		Bit_PackIPServiti_1=fopen("Bit_PackIPServiti_1.txt","w");
	}*/

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Inizializza_ATM(){

int i;

//printf("\n\n\n\t\t***** Inizializza_ATM *****\n"); // IL BUFFER IP E' IN BYTE

//MaxDimBuffer_ATM=ConvertiAdInteroSuperiore( ((MaxDimBuffer_IP[0]+MaxDimBuffer_IP[1])*UnByte/((double)(CellSize)*UnByte)) );  // IL BUFFER ATM E' IN CELLE ATM 
//MaxDimBuffer_ATM=	400	;  // IL BUFFER ATM E' IN CELLE ATM

DimBufferAttuale_ATM=0;
UltimoIstantePerdita_ATM=0.0;
Buffer_ATMVuoto=true;
InizioBusyPeriod_ATM=0.0;
ContributoIPA_ATM=0.0;
StimaIPADerivata_ATM=0.0;

for(i=0;i<NumBuffers_IP; i++){
	
	BusyPeriodAttivo[i]=false;
	UltimoIstantePerdita_IPoATM[i]=0.0;
	UltimoIstantePerdita_IPoATM_Delay[i]=0.0;
	InizioBusyPeriod_IPoATM[i]=0.0;
	StimaIPADerivata_IPoATM[i]=0.0;
	StimaIPADerivata_IPoATM_Delay[i]=0.0;
	ContributoIPA_IPoATM[i]=0.0;
	ContributoIPA_IPoATM_Delay[i]=0.0;
	PacchettiTotaliPersiNellaTramaATM_FlussoIP[i]=0.0;
	PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i]=0.0;
	TotSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]=0.0;
	TotQuadroSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]=0.0;
	TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[i]=0.0;

/////////////////////////////////////////////////////////// PID ///////////////////////////////////////////////////////////
	Kp[i]=3.0;
	Ki[i]=1.5;
	Kd[i]=1.25;
	e_k[i]=0.0;
	e_kMeno_1[i]=0.0;
	e_kMeno_2[i]=0.0;
	Ct_ATMSeparati[i]=Ct_ATM;
/////////////////////////////////////////////////////////// END PID ///////////////////////////////////////////////////////////
}

CelleArrivate=0.0;
CellePerse=0.0;
CelleServite=0.0;
BitArrivati_monitoraggioATM=0.0;
CampioneBitRate_monitoraggioATM=0.0;
CampioneBitRateQuadro_monitoraggioATM=0.0;

PacchettiTotaliPersiNellaTramaATM=0.0;

PPerdita_ATM=0.0;

InizializzaT_in_ATM();

if(GeneraConBurst){

	T_out_ATM =  (((double)(CellSize))*(double)(DimElementare)) / Ct_ATM;   //tempo per emettere sul canale ATM una cella ATM

}

/*printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Tempo per emettere sul canale ATM un byte= %e",( (double)(DimElementare) ) / Ct_ATM);
printf("\n\t****** Tempo che impiega la sorgente ATM[0] per emettere un byte= %e",T_in_ATM[0]/CellSize);
printf("\n\t****** Tempo che impiega la sorgente ATM[1] per emettere un byte= %e",T_in_ATM[1]/CellSize);
printf("\n\t****** Tempo per emettere sul canale ATM una cella= %e",T_out_ATM);
printf("\n\t****** Tempo sorgente ATM dal buffer IP 1 per emettere una cella= %e",T_in_ATM[0]);
printf("\n\t****** Tempo sorgente ATM dal buffer IP 2 per emettere una cella= %e",T_in_ATM[1]);
printf("\n");*/

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InizializzaT_in_ATM(){

int i;

/*if(BuffersInparallelo)
	T_in_ATM =  (((double)(CellSize))*(double)(DimElementare)) / Bp;
else*/
	for(i=0;i<NumBuffers_IP; i++)
		T_in_ATM[i] =  (((double)(CellSize))*(double)(DimElementare)) / Ct_IP[i];

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Inizializza_IP(){

int i;
track TracciaAttuale;

//printf("\n\n\n\t\t***** Inizializza_IP *****\n");

//Vincoli di SLA fissati
ReferenceIdealeLoss_0=1E-2; //2E-2 x il VoIP
ReferenceIdealeLoss_1=5E-3; 
//x la EqB
PPerdita_IP_IDEALE=ReferenceIdealeLoss_0;

MaxDelay[0]=20*1e-3;
MaxDelay[1]=20*1e-3;
ReferenceIdealeDelay_0=5E-2; 
ReferenceIdealeDelay_1=5E-2; 
//END Vincoli di SLA fissati

DimPacchettoMedio[0]=80*UnByte;
DimPacchettoMedio[1]=4534*UnByte;

for(i=0;i<NumBuffers_IP; i++){

DimBufferAttuale_IP[i]=0;

Buffer_IPVuoto[i]=true;

PPerdita_IP[i]=-1;

PacchettiArrivatiTot[i]=0.0;
PacchettiPersiTot[i]=0.0;
PacchettiServitiTot[i]=0.0;
Bit_PacchettiServitiTot[i]=0.0;
PacchettiServitioverATMTot[i]=0.0;

NumVolteSopraTarget[i]=0;
MediaDifferenzaSopraTarget[i]=0;
}

MaxDimBuffer_IP[0]=20*(DimPacchetto/UnByte);// IL BUFFER IP E' IN BYTE
MaxDimBuffer_IP[1]=10000*UnByte;// IL BUFFER IP E' IN BYTE

for(i=0; i<MaxNumConn_IP; i++){

        PacchettiArrivati[i]=0.0;
        PacchettiPersi[i]=0.0;
        /*DimMediaEffettivaDistrUniforme[i]=0.0;
        DimMediaEffettivaDistrEXP[i]=0.0;
		DimMediaEffettivaDistrBimodale[i]=0.0;
		PacchettiGrossiBimodale[i]=0;
		PacchettiPiccoliBimodale[i]=0;		

		a_PacchettiArrivatiDistrTrimodale[i]=0.0;
		b_PacchettiArrivatiDistrTrimodale[i]=0.0;
		c_PacchettiArrivatiDistrTrimodale[i]=0.0;*/

		PacchettiPersiNellaTramaATM[i]=0.0;
		PacchettiInRitardoNellaTramaATM[i]=0.0;
}

/////////////////////////////////////////////////// se genero secondo i pacchetti a dimensione variabile ///////////////////////////////////////////////////
if (!GeneraConBurst){

for(i=0;i<MaxNumConn_IP; i++)
        for(int j=0;j<=MaxDimensionePacchetto;j++)
                ContDimPackArrDistrUniforme[i][j]=0;

//double PerScaldarePareto=TempoInterarrivoProssimoPacchetto();//xche' il 1mo TempoInterarrivoProssimoPacchetto() secondo pareto e' un numero altissimo
//printf("\nStampo programmazione iniziale arrivo pacchetto per ciascuna delle %d sorgenti",(int)MaxNumConn_IP);
for(i=0; i<MaxNumConn_IP; i++){
        TracciaAttuale.evento = Arrivo_Pacchetto;
        TracciaAttuale.istante = 0.0;
        TracciaAttuale.Connessione = i;
        //printf("\nPer la connessione %d primo ArrivoPacchetto programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
}

L.sort();

//printf("\n\t****** DimensioneMediaPacchetto:\t%d byte",(int)DimensioneMediaPacchetto);
//printf("\n\t****** MaxDimensionePacchetto (caso Distr uniforme):\t%d byte",(int)MaxDimensionePacchetto);
/*printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps",Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps",Ct_IP[1]/((double)(UnMegabit)));
printf("\n\t****** Capacita' delle sorgenti IP: \t%lf Mbps",CtInSorgente_IP/((double)(UnMegabit)));
printf("\n\t****** Tempo per emettere sul canale IP un byte= %e",T_out_IP);
printf("\n\t****** Tempo che impiega una sorgente IP per emettere un byte= %e",T_in_IP);
printf("\n");*/

if(DistrUniformePack)
        printf("\nDistribuzione Uniforme nella dimensione del pacchetto");
if(DistrExpPack)
        printf("\nDistribuzione Esponenziale nella dimensione del pacchetto");
if(DistrBimodalePack){
        printf("\nDistribuzione Bimodale nella dimensione del pacchetto DimDistrBimodaleMIN=%d, DimDistrBimodaleMAX=%d",DimDistrBimodaleMIN,DimDistrBimodaleMAX);
		printf("\nProbabilita' di emettere un pacchetto di dimensione MIN=%lf",SogliaXDistrBimodale);
}
if(DistrTrimodalePack)
        printf("\nDistribuzione Trimodale nella dimensione del pacchetto");

}
/////////////////////////////////////////////////// END se genero secondo i pacchetti a dimensione variabile ///////////////////////////////////////////////////

/////////////////////////////////////////////////// se genero con burst ///////////////////////////////////////////////////
else{

MediaDurataBurst[0]= 1.0;			//VoIP	1.008, 1.587
//MediaDurataSilenzio[0]= 1.587;    
MediaDurataBurst[1]= 1.008;			//Video	10.0, 100.0
MediaDurataSilenzio[1]= 1.587;     

NumBurstGenerati=0.0;
DurataEffettivaDeiBurst=0.0;
NumSilenziGenerati=0.0;
DurataEffettivaDeiSilenzi=0.0;

NumConnAttive=0;

//Init_Bp[0]=16.0*UnKilobit;
Init_Bp[1]=0.001*UnKilobit;		//300.0

//MaxNumConn_IP_AttiveOra=MaxNumConn_IP	-	50;

for(i=0; i<MaxNumConn_IP; i++)
		Bp[i]=Init_Bp[CheBufferIPCorrente(i)];

for(i=0; i<MaxNumConn_IP; i++)
	T_in_IP[i]=( (double)(DimPacchetto) ) / Bp[i];

//MediaDurataBurstInCelle= MediaDurataBurst / T_in ;
////////////////////////////
if(GeneraSecondoClaudio){
        alfa=1.5;
		for(i=0;i<NumBuffers_IP; i++){
			deltaMediaDurataBurst[i]=MediaDurataBurst[i]*((alfa-1)/alfa);
			deltaMediaDurataSilenzio[i]=MediaDurataSilenzio[i]*((alfa-1)/alfa);
		}
}

////////////////////////////
//printf("\nStampo programmazione InizioBurst iniziali delle %d Connessioni",(int)MaxNumConn_IP_AttiveOra);

for(i=0; i<MaxNumConn_IP_AttiveOra; i++){
        TracciaAttuale.evento = Inizio_Burst;

        //TracciaAttuale.istante = 0.0;
        TracciaAttuale.istante = 0.0 + GenDurataSilenzio(i);

        //TracciaAttuale.istante = QuandoPrimoBurst[i];
        TracciaAttuale.Connessione = i;
        //printf("\nPer la connessione %d primo InizioBurst programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
		//printf("\nClock=%lf next burst per la conn %d al t=%lf",0.0,TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
}
/////////////////////////////////////////////////// END se genero con burst ///////////////////////////////////////////////////

//Generazione video
/*FileTracciatoVideo_Pacchetti=fopen("PacchettiVideo.txt","r");
FileTracciatoVideo_Tempi=fopen("TempiVideo.txt","r");
InizializzaTracciatoVideo();
fclose(FileTracciatoVideo_Pacchetti);
fclose(FileTracciatoVideo_Tempi);*/
//END Generazione video

//x file allocazione ideale
if(TecnicaRiallocazione==3 || TecnicaRiallocazione==5){

	Bit_PackIPServiti_0=fopen("Bit_PackIPServiti_0.txt","r");
	Bit_PackIPServiti_1=fopen("Bit_PackIPServiti_1.txt","r");							

	InizializzaTracciatoPacchettiIPServiti();

	fclose(Bit_PackIPServiti_0);
	fclose(Bit_PackIPServiti_1);
}
//END x file allocazione ideale

//x file allocazione IdealeLoss_2
if(TecnicaRiallocazione==5){

	PlossNellaTramaATM_0=fopen("PlossNellaTramaATM_0.txt","r");
	PlossNellaTramaATM_1=fopen("PlossNellaTramaATM_1.txt","r");							

	InizializzaTracciatoPloss();

	fclose(PlossNellaTramaATM_0);
	fclose(PlossNellaTramaATM_1);
}
//END x file allocazione IdealeLoss_2

L.sort();

/*Ct_IP[0]=(MaxNumConn_IP_AttiveOra-1) * (Init_Bp[0]*1.0);
Ct_IP[1]= 350.0* UnKilobit; //(MaxNumConn_IP_AttiveOra/2) * (Init_Bp[1]*0.28)	*/		

Ct_IP[0]=(MaxNumConn_IP_AttiveOra) * (Init_Bp[0]*1.0);
Ct_IP[1]= 1.0* UnKilobit;//2000.0* UnKilobit; 

//T_in_IP=( (double)(DimPacchetto) ) / Bp[0];
for(i=0;i<NumBuffers_IP; i++)
	T_out_IP[i]=( (double)(DimPacchetto) ) / Ct_IP[i];

/*printf("\n");
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps",Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps",Ct_IP[1]/((double)(UnMegabit)));
printf("\n\t****** Burstiness delle sorgenti IP[0]: \t%lf ",(1/(MediaDurataBurst[0]/(MediaDurataBurst[0]+MediaDurataSilenzio[0]))));
printf("\n\t****** Burstiness delle sorgenti IP[1]: \t%lf ",(1/(MediaDurataBurst[0]/(MediaDurataBurst[0]+MediaDurataSilenzio[1]))));
printf("\n\t****** Bp delle sorgenti IP[0]: \t%lf Mbps",Init_Bp[0]/((double)(UnMegabit)));
printf("\n\t****** Bp delle sorgenti IP[1]: \t%lf Mbps",Init_Bp[1]/((double)(UnMegabit)));
printf("\n\t****** Tempo per emettere dal buffer IP 1 una DimPacchetto= %e",T_out_IP[0]);
printf("\n\t****** Tempo per emettere dal buffer IP 2 una DimPacchetto= %e",T_out_IP[1]);
//printf("\n\t****** Tempo che impiega una sorgente IP per emettere una DimPacchetto= %e",T_in_IP);
printf("\n");*/

}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double DistribuzioneUniforme0_max(int max){

double NumeroRand0ToRandMax=rand();
double out=NumeroRand0ToRandMax/RAND_MAX;
	return out*max;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int GenAttualeDimPack(int NumConn){

int DimPacchettoAttuale;
double estrazione, DimdoublePacchettoAttuale;

if(DistrUniformePack){
        //estraggo DimdoublePacchettoAttuale che potra' essere x es: 0.19min 19.21MAX con MaxDimensionePacchetto=20
        DimdoublePacchettoAttuale=DistribuzioneUniforme0_max( (int)MaxDimensionePacchetto );
        //converto ad intero inferiore 0.19--->0 19.21--->19
        DimPacchettoAttuale = (int)DimdoublePacchettoAttuale;
        //aumento di 1 xche distribuzione uniforme del pacchetto tra 1 e MaxDimensionePacchetto; 0--->1(=minDimensionePacchetto) 19--->20(=MaxDimensionePacchetto)
        DimPacchettoAttuale++;
        if(DistrUniformePack)//senza questo if quando la distr exp genera una dim>MaxDimensionePacchetto la seguente istruzione da segmentation fault
                ContDimPackArrDistrUniforme[NumConn][DimPacchettoAttuale]++;
}

if(DistrExpPack){
		do
		{
			estrazione = NumeroRand(1.0);
		}
		while (estrazione==1.0);

	    DimdoublePacchettoAttuale=(double)DimensioneMediaPacchetto * -(log(1-estrazione));
        if(DimPacchettoAttuale==0)
                DimPacchettoAttuale++;
        DimMediaEffettivaDistrEXP[NumConn]+=DimPacchettoAttuale;
}

if(DistrBimodalePack){
	    estrazione=NumeroRand(1);
		if(estrazione <= SogliaXDistrBimodale){
			DimPacchettoAttuale=DimDistrBimodaleMIN;
			PacchettiPiccoliBimodale[NumConn]++;
		}
		else{
			DimPacchettoAttuale=DimDistrBimodaleMAX;
			PacchettiGrossiBimodale[NumConn]++;
		}
        DimMediaEffettivaDistrBimodale[NumConn]+=DimPacchettoAttuale;
}

if(DistrTrimodalePack){

	    estrazione=NumeroRand(1);    
		// Distribuzione Trimodale (a,b,c,pa,pb); papero Ajmone10/02 Trans. on Net. TRIMODAL(47 byte, 576 byte, 1500 byte; 0.559, 0.200)

		if(estrazione<=pa_DistrTrimodale){

			DimPacchettoAttuale=a_DimDistrTrimodale;
			a_PacchettiArrivatiDistrTrimodale[NumConn]++;
		}
		else{

			if ( (estrazione<=pa_DistrTrimodale + pb_DistrTrimodale) && (estrazione>pa_DistrTrimodale) ){
				
				DimPacchettoAttuale=b_DimDistrTrimodale;
				b_PacchettiArrivatiDistrTrimodale[NumConn]++;
			}
			else{
				
				DimPacchettoAttuale=c_DimDistrTrimodale;
				c_PacchettiArrivatiDistrTrimodale[NumConn]++;
			}
		}
}

return DimPacchettoAttuale;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void ArrivoPacchetto(int NumConn, double ClockAttuale){

track TracciaAttuale;
int DimPacchettoAttuale;

bool InviaCelle=1;

int BufferCorrente=CheBufferIPCorrente(NumConn);

/////////////////////////////////////////////////// se genero secondo i pacchetti a dimensione variabile ///////////////////////////////////////////////////
if (!GeneraConBurst){

	DimPacchettoAttuale = GenAttualeDimPack(NumConn);

	PacchettiArrivati[NumConn]++;

	if(DimBufferAttuale_IP[BufferCorrente]+DimPacchettoAttuale > MaxDimBuffer_IP[BufferCorrente]){
        PacchettiPersi[NumConn]++;
	}
	else{
        DimBufferAttuale_IP[BufferCorrente] += DimPacchettoAttuale;

        PacchettoInserito PacchettoInseritoAttuale;
        PacchettoInseritoAttuale.NumConn=NumConn;
        PacchettoInseritoAttuale.Dimensione=DimPacchettoAttuale;
        PacchettoInseritoAttuale.TInserimentoNelBuffer=ClockAttuale;
        CodaPacchettiInseriti[BufferCorrente].push(PacchettoInseritoAttuale);
        //printf("\nAl tempo %lf inserito per la conn %d un pack di dim %d DimBuffer=%d",ClockAttuale,PacchettoInseritoAttuale.NumConn, PacchettoInseritoAttuale.Dimensione,DimBufferAttuale_IP);
	}

	double TempoPerInserireIlCorrentePacchettoNelBuffer=T_in_IP[NumConn]*DimPacchettoAttuale;

	TracciaAttuale.evento = Arrivo_Pacchetto;
	TracciaAttuale.istante = ClockAttuale + TempoPerInserireIlCorrentePacchettoNelBuffer;
	TracciaAttuale.Connessione = NumConn;
	L.push_back(TracciaAttuale);

}
else{
	//
}
/////////////////////////////////////////////////// END se genero secondo i pacchetti a dimensione variabile ///////////////////////////////////////////////////

//Generazione video => se e' il buffer IP 1 e' il traffico video, altrimenti e' il VoIP
if(BufferCorrente==0){

	if(GeneraConBurst){// genera bursts, e' un if residuo che era in parallelo al if(!GeneraConBurst) di cui sopra

		// ******************** Genero con Burst il VoIP ********************
		if(ClockAttuale < IstanteFineAttualeBurst[NumConn]){

			DimPacchettoAttuale = (DimPacchetto/UnByte);
	
			PacchettiArrivati[NumConn]++;

			if(DimBufferAttuale_IP[BufferCorrente]+DimPacchettoAttuale > MaxDimBuffer_IP[BufferCorrente]){
				PacchettiPersi[NumConn]++;
			}
			else{
				DimBufferAttuale_IP[BufferCorrente] += DimPacchettoAttuale;

				PacchettoInserito PacchettoInseritoAttuale;
				PacchettoInseritoAttuale.NumConn=NumConn;
				PacchettoInseritoAttuale.Dimensione=DimPacchettoAttuale;
				PacchettoInseritoAttuale.TInserimentoNelBuffer=ClockAttuale;
				CodaPacchettiInseriti[BufferCorrente].push(PacchettoInseritoAttuale);
			}

			TracciaAttuale.evento = Arrivo_Pacchetto;
			TracciaAttuale.istante = ClockAttuale + T_in_IP[NumConn]; 
			TracciaAttuale.Connessione = NumConn;
			L.push_back(TracciaAttuale);
			//L.sort();
		}
		else { //if(ClockAttuale >= IstanteFineAttualeBurst[NumConn])

			TracciaAttuale.evento = Inizio_Burst;
			TracciaAttuale.istante = ClockAttuale + GenDurataSilenzio(NumConn);
        
			if (TracciaAttuale.istante==ClockAttuale)//puo' capitare con la distribuzione esponenziale nel GenDurataSilenzio e fa incartare il simulatore, 
				//( puo' darsi che sia anche il motivo dell'incarchiamento del simulatore accellerato )
				TracciaAttuale.istante+=T_in_IP[NumConn];
		
			TracciaAttuale.Connessione = NumConn;
			//printf("\nClock=%lf next burst per la conn %d al t=%lf",ClockAttuale,TracciaAttuale.Connessione,TracciaAttuale.istante);
			L.push_back(TracciaAttuale);
			//L.sort();

			NumConnAttive--;

			InviaCelle=false;
		}//END else //if(ClockAttuale > IstanteFineAttualeBurst[NumConn])
		
		// ******************** END Genero con Burst il VoIP ********************

	}//END if su (GeneraConBurst)
	else{
		//
	}

}//END if(BufferCorrente==0)
else{ //ossia BufferCorrente = 1, e' il traffico video

	// ******************** Genero video ********************
	double T_InterarrivoTracciatoVideo;
	
	DimPacchettoAttuale=(int)CampioniTracciatoVideo_Pacchetti[IndiceTracciatoVideo];
	T_InterarrivoTracciatoVideo=CampioniTracciatoVideo_Tempi[(IndiceTracciatoVideo+1)] * 1e-3;

	IndiceTracciatoVideo++;

	if(IndiceTracciatoVideo==(NumCampioniTracciatoVideo-1)){
		IndiceTracciatoVideo=0;
		T_InterarrivoTracciatoVideo=CampioniTracciatoVideo_Tempi[1];
	}
	
	PacchettiArrivati[NumConn]++;

	if(DimBufferAttuale_IP[BufferCorrente]+DimPacchettoAttuale > MaxDimBuffer_IP[BufferCorrente]){
        PacchettiPersi[NumConn]++;
	}
	else{
		DimBufferAttuale_IP[BufferCorrente] += DimPacchettoAttuale;

		PacchettoInserito PacchettoInseritoAttuale;
		PacchettoInseritoAttuale.NumConn=NumConn;
		PacchettoInseritoAttuale.Dimensione=DimPacchettoAttuale;
		PacchettoInseritoAttuale.TInserimentoNelBuffer=ClockAttuale;
		CodaPacchettiInseriti[BufferCorrente].push(PacchettoInseritoAttuale);
	}

	TracciaAttuale.evento = Arrivo_Pacchetto;
	TracciaAttuale.istante = ClockAttuale + T_InterarrivoTracciatoVideo; 
    TracciaAttuale.Connessione = NumConn;
    L.push_back(TracciaAttuale);
    //L.sort();
	// ******************** END Genero video ********************
}//END else su BufferCorrente 

if(Buffer_IPVuoto[BufferCorrente]){//se e' arrivato un pacchetto da qche parte risveglio l'InviaPacchettoNelCanale

	Buffer_IPVuoto[BufferCorrente]=false;

	//Risveglio l'invia pacchetto nel canale
	InviaPacchettoNelCanale(ClockAttuale, BufferCorrente);
}        

if(InviaCelle & BuffersInparallelo){

GeneraCelleDaPacchetto(ClockAttuale, DimPacchettoAttuale, NumConn);	

}//end del if (InviaCelle)

L.sort();

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void GeneraCelleDaPacchetto(double ClockAttuale, int DimPacchettoAttuale, int NumConn){

track TracciaAttuale;

double NumeroCelleATMNelPacchettoDouble=((DimPacchettoAttuale+4+12+8) / ((double)CellSize-(double)CellOverHead) ) ; 
int NumeroCelleATMNelPacchetto=ConvertiAdInteroSuperiore(NumeroCelleATMNelPacchettoDouble);

int Buffer_IP_Generante=CheBufferIPCorrente(NumConn);

Bit_PacchettiServitiTot[Buffer_IP_Generante]+=NumeroCelleATMNelPacchetto*(double)CellSize*(double)UnByte;

//printf("\n\n per il pacchetto arrivato al clock %lf programmo l'arrivo di %d celle ATM per i tempi",ClockAttuale,NumeroCelleATMNelPacchetto);

for( int k=0; k<NumeroCelleATMNelPacchetto; k++ ){

		TracciaAttuale.evento = Arrivo_Cella_ATM;
		TracciaAttuale.ConnIPGenerante=NumConn;

        TracciaAttuale.istante = ClockAttuale + k * T_in_ATM[Buffer_IP_Generante];

		if (k==0)
			TracciaAttuale.PrimaCellaDelPacchetto=true;
		else
			TracciaAttuale.PrimaCellaDelPacchetto=false;

		if (k==NumeroCelleATMNelPacchetto-1)
			TracciaAttuale.UltimaCellaDelPacchetto=true;
		else
			TracciaAttuale.UltimaCellaDelPacchetto=false;

		//printf("\n%lf",TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
}

// Qui non faccio la sort perche' la faccio comunque dopo aver chiamato tale funzione od in ArrivoPacchetto() od in InviiaPacchettoNelCanale()
//L.sort();
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int CheBufferIPCorrente(int NumConn){

int BufferCorrente;

if ( NumConn < (int)(MaxNumConn_IP_AttiveOra-1) )
	BufferCorrente=0;
else
	BufferCorrente=1;

//BufferCorrente=0;

//return BufferCorrente;
return 0;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InviaPacchettoNelCanale(double ClockAttuale, int Buffer_IP){

track TracciaAttuale;
double TempoEmmissioneAttualePacchetto;
int DimPacchettoAttualeDaEstrarre;
int ConnessioneIPInGioco;

if( DimBufferAttuale_IP[Buffer_IP] > 0 ){

        if(CodaPacchettiInseriti[Buffer_IP].empty())
                printf("\nATTENZIONE ERRUR in InviaPacchettoNelCanale!!! DimBufferAttuale_IP > 0 MA CodaPacchettiInseriti VUOTA");
        else{
                DimPacchettoAttualeDaEstrarre=CodaPacchettiInseriti[Buffer_IP].front().Dimensione;
				ConnessioneIPInGioco=CodaPacchettiInseriti[Buffer_IP].front().NumConn;

				if(GeneraConBurst)
					TempoEmmissioneAttualePacchetto=T_out_IP[Buffer_IP];
				else
					TempoEmmissioneAttualePacchetto=T_out_IP[Buffer_IP] * DimPacchettoAttualeDaEstrarre;
				
				PacchettiServitiTot[Buffer_IP]++;
				DimBufferAttuale_IP[Buffer_IP] -= DimPacchettoAttualeDaEstrarre;              
                //printf("\nAl tempo %lf servito per la conn %d un pack di dim %d inserito al tempo %lf DimBuffer=%d",ClockAttuale,CodaPacchettiInseriti.front().NumConn,DimPacchettoAttualeDaEstrarre,CodaPacchettiInseriti.front().TInserimentoNelBuffer,DimBufferAttuale_IP);
                //Ritardo_IP[CodaPacchettiInseriti[Buffer_IP].front().NumConn] += ClockAttuale + TempoEmmissioneAttualePacchetto -CodaPacchettiInseriti[Buffer_IP].front().TInserimentoNelBuffer;
                //Ritardo_IPTot += ClockAttuale + TempoEmmissioneAttualePacchetto - CodaPacchettiInseriti[Buffer_IP].front().TInserimentoNelBuffer;
                CodaPacchettiInseriti[Buffer_IP].pop();
        }
		
		if(!BuffersInparallelo && Buffer_IP==0){// ossia sono in cascata
			GeneraCelleDaPacchetto(ClockAttuale, DimPacchettoAttualeDaEstrarre, ConnessioneIPInGioco);
		}

		TracciaAttuale.evento = InviaPacchetto_NelCanale;
		TracciaAttuale.Buffer_IP = Buffer_IP;
		TracciaAttuale.istante = ClockAttuale + TempoEmmissioneAttualePacchetto;
		L.push_back(TracciaAttuale);
		L.sort();
}
else{
        
	//mi metto in attesa di un altro Busy Period per il buffer IP 
	//e, fin a quel momento, non rientrero' + in InviaPacchettoNelCanale() perche' non rischedulo il corrispondente evento
    Buffer_IPVuoto[Buffer_IP]=true;
}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void AggiornaContatoriTot_IP(){

int i;

for(i=0;i<NumBuffers_IP;i++){

	PacchettiArrivatiTot[i]=0.0;
	PacchettiPersiTot[i]=0.0;
}
	
for(i=0;i<MaxNumConn_IP;i++){

	int Buffer_IP=CheBufferIPCorrente(i);

    PacchettiPersiTot[Buffer_IP] += PacchettiPersi[i];
    PacchettiArrivatiTot[Buffer_IP] += PacchettiArrivati[i];
}
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void AggiornaContatoriTot_ATM(){

PacchettiTotaliPersiNellaTramaATM=0.0;

int i;

for(i=0;i<NumBuffers_IP;i++){
	PacchettiTotaliPersiNellaTramaATM_FlussoIP[i]=0.0;
	PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i]=0.0;
}

for(i=0;i<MaxNumConn_IP;i++){

	PacchettiTotaliPersiNellaTramaATM+=PacchettiPersiNellaTramaATM[i];
	
	int Buffer_IP=CheBufferIPCorrente(i);
    PacchettiTotaliPersiNellaTramaATM_FlussoIP[Buffer_IP] += PacchettiPersiNellaTramaATM[i];
	PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[Buffer_IP] += PacchettiInRitardoNellaTramaATM[i];
}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void CalcolaPPerdita_ATM(){

int i;

PPerdita_ATM = CellePerse / CelleArrivate;

double TotPacchettiArrivati=0.0;
for(i=0;i<NumBuffers_IP;i++)
	TotPacchettiArrivati+=PacchettiArrivatiTot[i];

PacketLossTramaATM=PacchettiTotaliPersiNellaTramaATM / TotPacchettiArrivati;

for(i=0;i<NumBuffers_IP;i++){
	PacketLossTramaATM_FlussoIP[i]=PacchettiTotaliPersiNellaTramaATM_FlussoIP[i] / PacchettiServitiTot[i];//PacchettiServitiTot con buff in cascata
	PacketDelayTramaATM_FlussoIP[i]=PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i] / PacchettiServitioverATMTot[i];
}

//Intervallo di confidenza loss

double zValidazPPerd=1.96;//1.96 oppure 2.58

double DevStandPopPPerd_ATM_FlussoIP[NumBuffers_IP];
double MaxScostamentoPPerd_ATM_FlussoIP[NumBuffers_IP];

for(i=0;i<NumBuffers_IP;i++){

	DevStandPopPPerd_ATM_FlussoIP[i]=( ( ((double)PacchettiTotaliPersiNellaTramaATM_FlussoIP[i]/(double)PacchettiServitiTot[i])*
	(1-((double)PacchettiTotaliPersiNellaTramaATM_FlussoIP[i]/(double)PacchettiServitiTot[i])) ) / (double)PacchettiServitiTot[i]);

	DevStandPopPPerd_ATM_FlussoIP[i]=sqrt(DevStandPopPPerd_ATM_FlussoIP[i]);
	MaxScostamentoPPerd_ATM_FlussoIP[i]=zValidazPPerd * DevStandPopPPerd_ATM_FlussoIP[i];

	if(PacketLossTramaATM_FlussoIP[i]>0.0){
		ScostamentoPercentualeLoss_IPoATM[i]= 
			( (PacketLossTramaATM_FlussoIP[i]+MaxScostamentoPPerd_ATM_FlussoIP[i])-PacketLossTramaATM_FlussoIP[i] ) / PacketLossTramaATM_FlussoIP[i];
		ScostamentoPercentualeLoss_IPoATM[i]*=100; 
	}
	else{
		ScostamentoPercentualeLoss_IPoATM[i]=0.0;
	}

}

//END Intervallo di confidenza loss

//Intervallo di confidenza Pdelay

double DevStandPopPdelay_ATM_FlussoIP[NumBuffers_IP];
double MaxScostamentoPdelay_ATM_FlussoIP[NumBuffers_IP];

for(i=0;i<NumBuffers_IP;i++){

	DevStandPopPdelay_ATM_FlussoIP[i]=( ( ((double)PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i]/(double)PacchettiServitioverATMTot[i])*
	(1-((double)PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i]/(double)PacchettiServitioverATMTot[i])) ) / (double)PacchettiServitioverATMTot[i]);

	DevStandPopPdelay_ATM_FlussoIP[i]=sqrt(DevStandPopPdelay_ATM_FlussoIP[i]);
	MaxScostamentoPdelay_ATM_FlussoIP[i]=zValidazPPerd * DevStandPopPdelay_ATM_FlussoIP[i];

	if(PacketDelayTramaATM_FlussoIP[i]>0.0){
		ScostamentoPercentualePdelay_IPoATM[i]= 
			( (PacketDelayTramaATM_FlussoIP[i]+MaxScostamentoPdelay_ATM_FlussoIP[i])-PacketDelayTramaATM_FlussoIP[i] ) / PacketDelayTramaATM_FlussoIP[i];
		ScostamentoPercentualePdelay_IPoATM[i]*=100; 
	}
	else{
		ScostamentoPercentualePdelay_IPoATM[i]=0.0;
	}

}

//END Intervallo di confidenza Pdelay

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void CalcolaPPerdita_IP(){

int i;

//per validazione statistica PPerdita_IP
/*double zValidazPPerd=2.58;//1.96 oppure 2.58
double gammaValidazPPerd=0.01;
double zDivgammaPPerd=0.0;
double zDivgammaPPerdQuadro=0.0;

double PerValidazPPerd;
//END per validazione statistica PPerdita_IP*/

for(i=0;i<NumBuffers_IP;i++)
	PPerdita_IP[i] = PacchettiPersiTot[i] / PacchettiArrivatiTot[i];

//per validazione statistica PPerdita_IP
/*zDivgammaPPerd=zValidazPPerd/gammaValidazPPerd;
zDivgammaPPerdQuadro=zDivgammaPPerd*zDivgammaPPerd;

PerValidazPPerd=( ((double)PacchettiArrivatiTot-(double)PacchettiPersiTot) / (double)PacchettiPersiTot);
DevStandPopPPerd=( ( ((double)PacchettiPersiTot/(double)PacchettiArrivatiTot)*(1-((double)PacchettiPersiTot/(double)PacchettiArrivatiTot)) ) / (double)PacchettiArrivatiTot);
DevStandPopPPerd=sqrt(DevStandPopPPerd);

MaxScostamentoPPerd=zValidazPPerd * DevStandPopPPerd;

if(PacchettiArrivatiTot > (PerValidazPPerd*zDivgammaPPerdQuadro) )
        StatisticaPPerdita_IPvalidata=true;
//END per validazione statistica PPerdita_IP*/
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void StampaRisultati(){

	StampaRisultati_IP();
	StampaRisultati_ATM();
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void StampaRisultati_IP(){

AggiornaContatoriTot_IP();
CalcolaPPerdita_IP();

//double MediaPPerdita_IPSingolaConn=0.0;

/*for(i=0;i<MaxNumConn_IP;i++){
        if(DistrUniformePack){
                for(int j=0;j<=MaxDimensionePacchetto;j++)
                        DimMediaEffettivaDistrUniforme[i] += ContDimPackArrDistrUniforme[i][j]*j;
                DimMediaEffettivaDistrUniforme[i] = DimMediaEffettivaDistrUniforme[i] / PacchettiArrivati[i];
        }
}*/

/*printf("\n********************************* Statistiche Connessione per Connessione *********************************");
for(i=0;i<MaxNumConn_IP;i++){
        printf("\n\n******** Per la connessione:\t%d ********",i);
        printf("\n\n numero pacchetti arrivati:\t%e",PacchettiArrivati[i]);
        printf("\n numero pacchetti persi :\t%e",PacchettiPersi[i]);
        printf("\n numero pacchetti serviti :\t%e",PacchettiServiti[i]);
        if(DistrExpPack)
                printf("\n DimMediaEffettivaDistrEXP:\t%lf byte",DimMediaEffettivaDistrEXP[i]/PacchettiArrivati[i]);
		if(DistrBimodalePack){
				printf("\n PacchettiPiccoliBimodale[%d]:\t%d",i,PacchettiPiccoliBimodale[i]);
				printf("\n PacchettiGrossiBimodale[%d]:\t%d",i,PacchettiGrossiBimodale[i]);
                printf("\n DimMediaEffettivaDistrBimodale:\t%lf byte",DimMediaEffettivaDistrBimodale[i]/PacchettiArrivati[i]);
		}
        if(DistrUniformePack){
                printf("\n DimMediaEffettivaDistrUniforme:\t%lf byte",DimMediaEffettivaDistrUniforme[i]);
                printf("\nElenco dei pacchetti arrivati");
                for(int j=0;j<=MaxDimensionePacchetto;j++)
                        printf("\npacchetti arrivati di dimensione %d byte:\t%d",j,ContDimPackArrDistrUniforme[i][j]);
        }
		if(DistrTrimodalePack){

			double pacchettitotaliarrivati=a_PacchettiArrivatiDistrTrimodale[i]+b_PacchettiArrivatiDistrTrimodale[i]+c_PacchettiArrivatiDistrTrimodale[i];
			printf("\n\npacchetti totali arrivati Conn%d=%lf",i,pacchettitotaliarrivati);
			printf("\na_PacchettiArrivatiDistrTrimodale=%e",a_PacchettiArrivatiDistrTrimodale[i]);
			printf("\nb_PacchettiArrivatiDistrTrimodale=%e",b_PacchettiArrivatiDistrTrimodale[i]);
			printf("\nc_PacchettiArrivatiDistrTrimodale=%e",c_PacchettiArrivatiDistrTrimodale[i]);
			printf("\npa_DistrTrimodale=%lf",a_PacchettiArrivatiDistrTrimodale[i]/pacchettitotaliarrivati);
			printf("\npb_DistrTrimodale=%lf",b_PacchettiArrivatiDistrTrimodale[i]/pacchettitotaliarrivati);
			printf("\npc_DistrTrimodale=%lf",c_PacchettiArrivatiDistrTrimodale[i]/pacchettitotaliarrivati);	
		}
}
printf("\n\n********************************* END Statistiche Connessione per Connessione *********************************");*/

if (GeneraConBurst){

	printf("\n\nNumConnAttive=%d",NumConnAttive);
	
	printf("\n\nNumBurstGenerati=%d",(int)NumBurstGenerati);
	printf("\nDurataEffettivaDeiBurst =%lf",(DurataEffettivaDeiBurst/NumBurstGenerati));

	printf("\n\nNumSilenziGenerati=%d",(int)NumSilenziGenerati);
	printf("\nDurataEffettivaDeiSilenzi =%lf",(DurataEffettivaDeiSilenzi/NumSilenziGenerati));
}
	
	
printf("\n\nnumero pacchetti arrivati in totale %e",PacchettiArrivatiTot);
printf("\nnumero pacchetti persi in totale %e",PacchettiPersiTot);

printf("\n\n\t****** PPerdita_IP=%e ******",PPerdita_IP);
//printf("\n\t****** Ritardo_IP Medio su tutte le celle %lf ******",Ritardo_IPTot/PacchettiServitiTot);
//printf("\n\t****** DimBufferAttuale_IP=%d byte******",DimBufferAttuale_IP);

//printf("\n\n\t****** Bit Arrivati=%e ******",bitArrivati_IP);
/*printf("\n\t****** Bit Persi=%e ******",bitPersi_IP);
printf("\n\t****** PPerdita_IPInBit=%e ******",(bitPersi_IP/bitArrivati_IP));*/

//for(i=0;i<MaxNumConn_IP;i++)
//        printf("\n\tRitardo_IP Medio su tutte le celle della connessione %d:\t%lf",i,Ritardo_IP[i]/PacchettiServiti[i]);

/*for(i=0;i<MaxNumConn_IP;i++){
        PPerdita_IPSingolaConn[i]= PacchettiPersi[i] / PacchettiArrivati[i];
        MediaPPerdita_IPSingolaConn += PPerdita_IPSingolaConn[i];
        //printf("\nPPerdita_IP connessione %d =\t%lf",i,PPerdita_IPSingolaConn[i]);
        DiffRispettoPPerdita_IP[i] = PPerdita_IPSingolaConn[i] - PPerdita_IP;
}

MediaPPerdita_IPSingolaConn=MediaPPerdita_IPSingolaConn/((double)MaxNumConn_IP);*/

//printf("\n\nMediaPPerdita_IPSingolaConn=\t%lf\n",MediaPPerdita_IPSingolaConn);

//for(i=0;i<MaxNumConn_IP;i++)
//        printf("\nDiffRispettoPPerdita_IP connessione %d =\t%lf",i,(DiffRispettoPPerdita_IP[i]));

//printf("\n\n\tLa minima DiffRispettoPPerdita_IP =\t%lf",(TrovaMinimo(DiffRispettoPPerdita_IP,MaxNumConn_IP)));
//printf("\n\n\tLa Massima DiffRispettoPPerdita_IP =\t%lf",(TrovaMassimo(DiffRispettoPPerdita_IP,MaxNumConn_IP)));

/*if(StatisticaPPerdita_IPvalidata)
        printf("\n\n***** LA STATISTICA SULLA MISURA DELLA PPerdita_IP E VALIDATA *****");
else
        printf("\n\n***** LA STATISTICA SULLA MISURA DELLA PPerdita_IP NON E VALIDATA *****");*/

/*printf("\nDeviazione standard Popolazione PPerd e' %lf",DevStandPopPPerd);
printf("\nScostamento teorico PPerd e' %lf",MaxScostamentoPPerd);
printf("\nMAX Scostamento sopra PPerd previsto teoricamente:%lf",PPerdita_IP+MaxScostamentoPPerd);
printf("\nMAX Scostamento sotto PPerd previsto teoricamente:%lf",PPerdita_IP-MaxScostamentoPPerd);
printf("\nOcio il valore VERO PPerd e' %lf + o - %lf",PPerdita_IP,MaxScostamentoPPerd);*/

/*if(DistrUniformePack)
        printf("\n\n.....e la Distribuzione Uniforme nella dimensione del pacchetto");
if(DistrExpPack)
        printf("\n\n.....e la Distribuzione Esponenziale nella dimensione del pacchetto");
if(DistrBimodalePack)
        printf("\n\n.....e la Distribuzione Bimodale nella dimensione del pacchetto");
if(DistrTrimodalePack)
        printf("\n\n.....e la Distribuzione Trimodale nella dimensione del pacchetto");*/

printf("\n");
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void StampaRisultati_ATM(){


AggiornaContatoriTot_ATM();
CalcolaPPerdita_ATM();

printf("\n\nnumero celle arrivate in totale %e", CelleArrivate);
printf("\nnumero celle perse in totale %e", CellePerse);
printf("\nnumero celle servite in totale %e", CelleServite);

//ATTENZIONE il numero di bit ATM arrivati in totale a meno dell'Overhead ATM sara' un po' piu' grande dei bit IP arrivati in quanto faccio la conversione:
//int NumeroCelleATMNelPacchetto=( (int) (DimPacchettoAttuale/ (CellSize-CellOverHead) ) ) + 1.0; ...
//...che spreca un po' di bit nella trama ATM
//printf("\n\nnumero bit ATM arrivati in totale a meno dell'Overhead ATM %e", CelleArrivate*(CellSize-CellOverHead)*DimElementare);
/*printf("\nnumero bit ATM persi in totale compreso l'OverHead ATM %e", CellePerse*CellSize*DimElementare);
printf("\nnumero bit ATM serviti in totale compreso l'OverHead ATM %e", CelleServite*CellSize*DimElementare);*/

printf("\n\n\t****** PPerdita_ATM=%e ******", PPerdita_ATM);
//printf("\n\t****** PPerdita_ATMInBit=%e ******",((CellePerse*CellSize*DimElementare)/(CelleArrivate*CellSize*DimElementare)));

/*for (i=0;i<MaxNumConn_IP;i++){

	printf("\n\nPacchettiPersiNellaTramaATM[%d]=%e",i,PacchettiPersiNellaTramaATM[i]);
	printf("\nPacketLossTramaATM[%d]=%e",i,(PacchettiPersiNellaTramaATM[i]/PacchettiArrivati[i]));
}*/

printf("\n\nPacketLossTramaATM =%e",PacketLossTramaATM);
printf("\n\nPacketLossTramaATM_FlussoIP[0] =%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1] =%e",PacketLossTramaATM_FlussoIP[1]);

printf("\n\nPacketDelayTramaATM_FlussoIP[0] =%e",PacketDelayTramaATM_FlussoIP[0]);
printf("\nPacketDelayTramaATM_FlussoIP[1] =%e",PacketDelayTramaATM_FlussoIP[1]);

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void StampaLoss(){
	
	printf("\n\n\n");

	printf("\nPacketLossTramaATM=%e",PacketLossTramaATM);
	printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
	printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);
	printf("\nPPerdita_IP[%d]=%e",0,PPerdita_IP[0]);
	printf("\nPPerdita_IP[%d]=%e",1,PPerdita_IP[1]);

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Svuota_EventiRimastiInCoda(){

//printf("\n\t***** Ora stampo gli eventi rimasti nella coda degli eventi *****");
while (!L.empty()){
	L.pop_front();
}//fine stampa eventi rimasti nella coda degli eventi */
//printf("\n");
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int TrovaIndiceMinimo(double *Array, int MaxDim){

	int res;
	
	double minimo=1000000000000.0;
	for(int i=0;i<MaxDim;i++)
		if(Array[i]<minimo){
			minimo=Array[i];
			res=i;
		}

	return res;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double TrovaMinimo(double *Array, int MaxDim){

double minimo=1000000000000.0;
for(int i=0;i<MaxDim;i++)
        if(Array[i]<minimo)
                minimo=Array[i];
return minimo;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double TrovaMassimo(double *Array, int MaxDim){

double Max=-1000000000000.0;
for(int i=0;i<MaxDim;i++)
        if(Array[i]>Max)
                Max=Array[i];
return Max;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void AzzeraContatoriPerRiprendereIlContoDalRegime(){

	//printf("\n\n\t\t***** AzzeraContatoriPerRiprendereIlContoDalRegime *****\n");

	// 20.03.12 - OCIO: mi ero dimenticato il seguente azzeramente che mi faceva sballare la media calcolata in addestramento rispetto a quella nel controllo
	BitArrivati_monitoraggioATM=0.0;

	DimBufferAttuale_ATM = DimBufferAttuale_ATM_scelto ;

int i;

////////////////// Lato IP //////////////////

for(i=0;i<NumBuffers_IP; i++){
	PacchettiArrivatiTot[i]=0.0;
	PacchettiPersiTot[i]=0.0;
	PacchettiServitiTot[i]=0.0;
	Bit_PacchettiServitiTot[i]=0.0;
	PacchettiServitioverATMTot[i]=0.0;
}

for(i=0; i<MaxNumConn_IP; i++){

        PacchettiArrivati[i]=0.0;
        PacchettiPersi[i]=0.0;        
        //DimMediaEffettivaDistrUniforme[i]=0.0;
        //DimMediaEffettivaDistrEXP[i]=0.0;
		//DimMediaEffettivaDistrBimodale[i]=0.0;
		//PacchettiGrossiBimodale[i]=0;
		//PacchettiPiccoliBimodale[i]=0;		

		//a_PacchettiArrivatiDistrTrimodale[i]=0.0;
		//b_PacchettiArrivatiDistrTrimodale[i]=0.0;
		//c_PacchettiArrivatiDistrTrimodale[i]=0.0;

		PacchettiPersiNellaTramaATM[i]=0.0;
		PacchettiInRitardoNellaTramaATM[i]=0.0;
}

//for(i=0;i<MaxNumConn_IP; i++)
//        for(int j=0;j<=MaxDimensionePacchetto;j++)
//                ContDimPackArrDistrUniforme[i][j]=0;

////////////////// END Lato IP //////////////////

////////////////// Lato ATM //////////////////

CelleArrivate=0.0;
CellePerse=0.0;
CelleServite=0.0;

for(i=0;i<NumBuffers_IP;i++){
	PacchettiTotaliPersiNellaTramaATM_FlussoIP[i]=0.0;
	PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[i]=0.0;
}

PacchettiTotaliPersiNellaTramaATM=0.0;

////////////////// END Lato ATM //////////////////
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
//funzione di ordinamento della lista eventi 
bool track::operator < (track T){

	return (istante<T.istante);
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double NumeroRand(int max){

  double tmpnum;

  double init = rand();

  tmpnum = init / RAND_MAX;
  tmpnum = tmpnum * max;

  return tmpnum;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void ArrivoCella_ATM(double ClockAttuale, int ConnIPGenerante, bool PrimaCellaDelPacchetto, bool UltimaCellaDelPacchetto){

int i;
	
if(PrimaCellaDelPacchetto)
	PackInCorsoTramaATMGiaPerso[ConnIPGenerante]=false;

int FlussoIPCorrente=CheBufferIPCorrente(ConnIPGenerante);

CelleArrivate++;
BitArrivati_monitoraggioATM+=(double)CellSize*(double)UnByte;

if(DimBufferAttuale_ATM+1 > MaxDimBuffer_ATM){

	CellePerse++;
	
	//IPA relativa al flusso ATM in se per se
	//UltimoIstantePerdita_ATM=ClockAttuale;

	//IPA relativa al flusso IPoATM
	if (!PackInCorsoTramaATMGiaPerso[ConnIPGenerante]){

		//la prima cella persa che implica l'ultimo pacchetto perso definisce UltimoIstantePerdita_IPoATM
		//le altre celle ATM che si perderanno devono avere effetto solo sulla UltimoIstantePerdita_ATM
		UltimoIstantePerdita_IPoATM[FlussoIPCorrente]=ClockAttuale;
		
		PacchettiPersiNellaTramaATM[ConnIPGenerante]++;
		PackInCorsoTramaATMGiaPerso[ConnIPGenerante]=true;
	}
}
else{
	DimBufferAttuale_ATM++;

	//Calcolo PDelay
	if(UltimaCellaDelPacchetto && PackInCorsoTramaATMGiaPerso[ConnIPGenerante]==false)//se il pack è entrato tutto ne calcolo il delay
	{
		
		PacchettiServitioverATMTot[FlussoIPCorrente]++;
		
		DelayPackAttuale=(DimBufferAttuale_ATM*(double)CellSize*(double)UnByte)/Ct_ATM;
		/*if(FlussoIPCorrente==0){
			printf("\nDelayPackAttuale=%lf, %lf[ms]\t\t\t", DelayPackAttuale, DelayPackAttuale*1e3); 
			system("pause");
		}*/
		UltimoDelayFlussoIPoATM[FlussoIPCorrente]=DelayPackAttuale;

		if(DelayPackAttuale>=MaxDelay[FlussoIPCorrente]){
			
			if(TecnicaRiallocazione!=6){
				UltimoIstantePerdita_IPoATM_Delay[FlussoIPCorrente]=ClockAttuale;
				PacchettiInRitardoNellaTramaATM[ConnIPGenerante]++;
			}
			else{//if(TecnicaRiallocazione!=6)
				
				if(NumeroRand(1)>ReferenceIdealeDelay_0){
					//rialloco la banda x ottenere l'esatto delay previsto;
					Ct_ATM=(DimBufferAttuale_ATM*(double)CellSize*(double)UnByte)/MaxDelay[FlussoIPCorrente];

					EntrateInRiallocazioneIdealDelay++;

					Ct_ATM_Tot+=Ct_ATM;
					Ct_ATM_TotQuadro+=Ct_ATM*Ct_ATM;
				}
				else{
					UltimoIstantePerdita_IPoATM_Delay[FlussoIPCorrente]=ClockAttuale;
					PacchettiInRitardoNellaTramaATM[ConnIPGenerante]++;
				}
			
			}//END else if(TecnicaRiallocazione!=6)
		
		}//END if(DelayPackAttuale>=MaxDelay[FlussoIPCorrente])
	}
	//Calcolo PDelay

	/*CellaInserita CellaInseritaAttuale;
	CellaInseritaAttuale.NumConn=ConnIPGenerante;
	CellaInseritaAttuale.TInserimentoNelBuffer=ClockAttuale;
	CodaCelleInserite.push(CellaInseritaAttuale);*/
}

if(Buffer_ATMVuoto){

	Buffer_ATMVuoto=false;

	BusyPeriodAttivo[FlussoIPCorrente]=true; 

	//Comincia un nuovo Busy Period

	//IPA relativa al flusso ATM in se per se
	//InizioBusyPeriod_ATM=ClockAttuale; 

	//IPA relativa al flusso IPoATM
	InizioBusyPeriod_IPoATM[FlussoIPCorrente]=ClockAttuale;

	//Risveglio l'invia cella nel canale
	InviaCellaNelCanale_ATM(ClockAttuale);          
}
else{ //se il buffer non e' vuoto...

	// ... ma non e' ancora arrivato un pacchetto di quel specifico flusso IP...
	// ... aggiorno il suo InizioBusyPeriod_IPoATM
	if(BusyPeriodAttivo[FlussoIPCorrente]==false){

		InizioBusyPeriod_IPoATM[FlussoIPCorrente]=ClockAttuale;
		BusyPeriodAttivo[FlussoIPCorrente]=true; 
	}

}// end else if(Buffer_ATMVuoto)
    
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InviaCellaNelCanale_ATM(double ClockAttuale){

track TracciaAttuale;

int i;

if(DimBufferAttuale_ATM > 0){
        
	//ConnessioneIPInGioco=CodaCelleInserite.front().NumConn;
    //CodaCelleInserite.pop();
	
	DimBufferAttuale_ATM--;
    CelleServite++;

	TracciaAttuale.evento = InviaCella_NelCanale_ATM;
	TracciaAttuale.istante = ClockAttuale + T_out_ATM;
	L.push_back(TracciaAttuale);
	L.sort();
}
else{//il buffer si e' svuotato

	//fine dell'ultimo Busy Period
	//aggiorno ContributoIPA se ho perso durante il Busy Period
	
	//questo riguarda l'IPA slegata dall'IPoATM
	//if(UltimoIstantePerdita_ATM!=0.0)
		//ContributoIPA_ATM += ( UltimoIstantePerdita_ATM - InizioBusyPeriod_ATM );
	
	//questo riguarda l'IPA dell'IPoATM
	for(i=0;i<NumBuffers_IP; i++){
		if(UltimoIstantePerdita_IPoATM[i]!=0.0)
			ContributoIPA_IPoATM[i] += ( UltimoIstantePerdita_IPoATM[i] - InizioBusyPeriod_IPoATM[i] );
		if(UltimoIstantePerdita_IPoATM_Delay[i]!=0.0)
			ContributoIPA_IPoATM_Delay[i] += ( UltimoIstantePerdita_IPoATM_Delay[i] - InizioBusyPeriod_IPoATM[i] );
	}

	//riinizializzo
	//InizioBusyPeriod_ATM=0.0; 
	//UltimoIstantePerdita_ATM=0.0;
	
	for(i=0;i<NumBuffers_IP; i++){

		BusyPeriodAttivo[i]=false; 
		InizioBusyPeriod_IPoATM[i]=0.0; 
		UltimoIstantePerdita_IPoATM[i]=0.0;
		UltimoIstantePerdita_IPoATM_Delay[i]=0.0;
	}
        
	//mi metto in attesa di un altro Busy Period 
	//e, fin a quel momento, non rinetrero' + in InviaCellaNelCanale() perche' non rischedulo il corrispondente evento
    Buffer_ATMVuoto=true;
	//printf("\n\n\t***** ClokAttuale=%lf Messo a sopire l'InviaCellaNelCanale per la stazione %d",ClockAttuale,Stazione);
}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM(double ClockAttuale){

	int i;

AggiornaContatoriTot_IP();
AggiornaContatoriTot_ATM();
CalcolaPPerdita_IP();
CalcolaPPerdita_ATM();

for(i=0;i<NumBuffers_IP; i++)
	ScartoQuadraticoMedio_IPTramaATM += pow( ( PacchettiPersiTot[i]-PacchettiTotaliPersiNellaTramaATM_FlussoIP[i] ), 2);

if (BuffersInparallelo){
	//
}
else{	
	
	LossRate_IDEALI_0=PacchettiServitiTot[0]*DimPacchettoMedio[0]*ReferenceIdealeLoss_0/DimIntervalloRiallocazioniBanda;//in bps
	LossRate_IDEALI_1=PacchettiServitiTot[1]*DimPacchettoMedio[1]*ReferenceIdealeLoss_1/DimIntervalloRiallocazioniBanda;

	DelayRate_IDEALI_0=PacchettiServitioverATMTot[0]*DimPacchettoMedio[0]*ReferenceIdealeDelay_0/DimIntervalloRiallocazioniBanda;//in bps
	DelayRate_IDEALI_1=PacchettiServitioverATMTot[1]*DimPacchettoMedio[1]*ReferenceIdealeDelay_1/DimIntervalloRiallocazioniBanda;
}

////////////////////////////////////////////////
//scelta tecnica di riallocazione
if(TecnicaRiallocazione==0)
	RiallocazioneCt_ATM_IPA(ClockAttuale);
if(TecnicaRiallocazione==1)
	RiallocazioneCt_ATM_BEqBufferless(ClockAttuale);
if(TecnicaRiallocazione==2)
	RiallocazioneCt_ATM_PID(ClockAttuale);
if(TecnicaRiallocazione==3)
	RiallocazioneCt_ATM_Ideale(ClockAttuale);
if(TecnicaRiallocazione==4)
	RiallocazioneCt_ATM_Bm(ClockAttuale);
if(TecnicaRiallocazione==5)
	RiallocazioneCt_ATM_IdealeLoss_2(ClockAttuale);
if(TecnicaRiallocazione==6)
	RiallocazioneCt_ATM_Bm(ClockAttuale);
//END scelta tecnica di riallocazione
////////////////////////////////////////////////

//contazeri per RCBC
/*
if(TecnicaRiallocazione==0){

	if(PacketLossTramaATM_FlussoIP[1] == 0.0)
		contazeri++;
	//EURASIP: if zero values of PDelay are registered for six consecutive times, the bandwidth is decreased of 2%
	if(contazeri == 6){
		Ct_ATM = Ct_ATM - (Ct_ATM*0.02);
		contazeri=0;
	}

}*/
//END contazeri per RCBC

/////////// Stampa situazione ///////////

printf("\n");
for(i=0;i<NumBuffers_IP; i++)
	printf("\nPPerdita_IP[%d]=%e",i,PPerdita_IP[i]);
printf("\n");
//printf("\nPacketLossTramaATM=%e",PacketLossTramaATM);
printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);
/* esempio di delay
buffer 150400 bit (100 celle DVB)
Ct_ATM 829000bps, .829 kbps
0.181s 181ms di ritardo
 END esempio di delay */
printf("\nPacketDelayTramaATM_FlussoIP[0]=%e",PacketDelayTramaATM_FlussoIP[0]);
printf("\nPacketDelayTramaATM_FlussoIP[1]=%e",PacketDelayTramaATM_FlussoIP[1]);

/* e' relativo allo scostamento [%] previsto rispetto all'intervallo di confidenza:
printf("\n\nScostamentoPercentualeLoss_IPoATM[0]=%lf",ScostamentoPercentualeLoss_IPoATM[0]);
printf("\nScostamentoPercentualeLoss_IPoATM[1]=%lf",ScostamentoPercentualeLoss_IPoATM[1]);
printf("\nScostamentoPercentualePdelay_IPoATM[0]=%lf",ScostamentoPercentualePdelay_IPoATM[0]);
printf("\nScostamentoPercentualePdelay_IPoATM[1]=%lf",ScostamentoPercentualePdelay_IPoATM[1]);*/

if(TecnicaRiallocazione==0){
	//printf("\nSpintaGradiente=%lf Kbps", SpintaGradiente/((double)(UnKilobit)));
	printf("\n\nStimaIPADerivata_IPoATM[0]=%lf Kbps", StimaIPADerivata_IPoATM[0]/((double)(UnKilobit)));
	printf("\nStimaIPADerivata_IPoATM[1]=%lf Kbps", StimaIPADerivata_IPoATM[1]/((double)(UnKilobit)));
	printf("\n\nStimaIPADerivata_IPoATM_Delay[0]=%lf Kbps", StimaIPADerivata_IPoATM_Delay[0]/((double)(UnKilobit)));
	printf("\nStimaIPADerivata_IPoATM_Delay[1]=%lf Kbps", StimaIPADerivata_IPoATM_Delay[1]/((double)(UnKilobit)));

	printf("\n\nContributoIPA_IPoATM[0]=%lf ", ContributoIPA_IPoATM[0]);
	printf("\nContributoIPA_IPoATM[1]=%lf ", ContributoIPA_IPoATM[1]);
	printf("\nContributoIPA_IPoATM_Delay[0]=%lf ", ContributoIPA_IPoATM_Delay[0]);
	printf("\nContributoIPA_IPoATM_Delay[1]=%lf ", ContributoIPA_IPoATM_Delay[1]);
}

printf("\n\nLossRateNellaTramaATM_FlussoIP[0]=%lf Kbps", (PacchettiTotaliPersiNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nLossRateNellaTramaATM_FlussoIP[1]=%lf Kbps", (PacchettiTotaliPersiNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nLossRate_IDEALE_0=%lf Kbps", LossRate_IDEALI_0/((double)(UnKilobit)));
printf("\nLossRate_IDEALE_1=%lf Kbps", LossRate_IDEALI_1/((double)(UnKilobit)));
printf("\n");
printf("\nDelayRateNellaTramaATM_FlussoIP[0]=%lf Kbps", (PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nDelayRateNellaTramaATM_FlussoIP[1]=%lf Kbps", (PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nDelayRate_IDEALE_0=%lf Kbps", DelayRate_IDEALI_0/((double)(UnKilobit)));
printf("\nDelayRate_IDEALE_1=%lf Kbps", DelayRate_IDEALI_1/((double)(UnKilobit)));
printf("\n");

//StampaRisultati();

//Ct_ATM+=ProvaIncremetoCt_ATM;
T_out_ATM =  ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;
printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps", Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps", Ct_IP[1]/((double)(UnMegabit)));
/*printf("\n\t****** Simulated CellTax: \t%lf ", ( ((Ct_ATM-(Ct_IP[0]+Ct_IP[1]))/(Ct_IP[0]+Ct_IP[1])) * 100) );
printf("\n\n\t****** CellTax teorico: \t%lf ", CellTaxTeorico);
printf("\n\t****** Capacita' del canale ATM secondo CellTax teorico: \t%lf Mbps", Ct_ATMSecondoCellTaxTeorico/((double)(UnMegabit)));*/
//printf("\n\t****** Tempo per emettere sul canale ATM una cella= %e",T_out_ATM);

printf("\n\n"); 
printf("\nUltimoDelayFlussoIPoATM[0]=%lf, %lf[ms]\t\t\t", UltimoDelayFlussoIPoATM[0], UltimoDelayFlussoIPoATM[0]*1e3); 
printf("\nUltimoDelayFlussoIPoATM[1]=%lf, %lf[ms]\t\t\t", UltimoDelayFlussoIPoATM[1], UltimoDelayFlussoIPoATM[1]*1e3); 

//system("pause");

/////////// END Stampa situazione ///////////

////////////////////////////////////////////////
////////////////////////////////////////////////
//calcoli relativi alle misure medie su tutta la simulazione dell'intervallo di confidenze etc.

EntrateInRiallocazione++;
if(TecnicaRiallocazione==6)
	EntrateInRiallocazioneIdealDelay++;

Ct_ATM_Tot+=Ct_ATM;
Ct_ATM_TotQuadro+=Ct_ATM*Ct_ATM;

if(!LavoraSuPDelay){//uso un unico contatore totale sulle medie su tutta la simulazione col nome riferito alla loss, ma sulla base di questo booleano ci posso mettere anche il delay

	for(i=0;i<NumBuffers_IP; i++){

		if(PacketLossTramaATM_FlussoIP[i]>0.0){//x evitare stampe di infinito
			TotSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]+=PacketLossTramaATM_FlussoIP[i];
			TotQuadroSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]+=PacketLossTramaATM_FlussoIP[i]*PacketLossTramaATM_FlussoIP[i];
		}

		TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[i]+=ScostamentoPercentualeLoss_IPoATM[i];		
	}

	if(PacketLossTramaATM_FlussoIP[0]>ReferenceIdealeLoss_0){
		NumVolteSopraTarget[0]++;
		MediaDifferenzaSopraTarget[0]+=(PacketLossTramaATM_FlussoIP[0]-ReferenceIdealeLoss_0);
	}
	if(PacketLossTramaATM_FlussoIP[1]>ReferenceIdealeLoss_1){
		NumVolteSopraTarget[1]++;
		MediaDifferenzaSopraTarget[1]+=(PacketLossTramaATM_FlussoIP[1]-ReferenceIdealeLoss_1);
	}

}//END if(!LavoraSuPDelay)
else
{
	
	for(i=0;i<NumBuffers_IP; i++){

		if(PacketDelayTramaATM_FlussoIP[i]>0.0){//x evitare stampe di infinito
			TotSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]+=PacketDelayTramaATM_FlussoIP[i];
			TotQuadroSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]+=PacketDelayTramaATM_FlussoIP[i]*PacketDelayTramaATM_FlussoIP[i];
		}

		TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[i]+=ScostamentoPercentualePdelay_IPoATM[i];		
	}

	if(PacketDelayTramaATM_FlussoIP[0]>ReferenceIdealeDelay_0){
		NumVolteSopraTarget[0]++;
		MediaDifferenzaSopraTarget[0]+=(PacketDelayTramaATM_FlussoIP[0]-ReferenceIdealeDelay_0);
	}
	if(PacketDelayTramaATM_FlussoIP[1]>ReferenceIdealeDelay_1){
		NumVolteSopraTarget[1]++;
		MediaDifferenzaSopraTarget[1]+=(PacketDelayTramaATM_FlussoIP[1]-ReferenceIdealeDelay_1);
	}
}

/*
//non appena mi è andato sotto il target considero la convergenza raggiunta e ricomincio a calcolare la media e varianza delle misure sull'intera simulazione
//OCIO!: se mi ha raggiunto il target, ma poi oscilla sopra e sotto eccessivamente me ne accorgerò dalla media e varianza sull'intera simulazione dopo la convergenza
if(PacketLossTramaATM_FlussoIP[0]<ReferenceIdealeLoss_0){

	//lo disabilito se non voglio azzerare i contatori dopo l'arrivo in convergenza (cosa che facevo x la subsection I dei risultati, quella che visualizza l'intervallo di confidenza)
	//AggiornamentoPerArrivoInConvergenza=true;

}

if(AggiornamentoPerArrivoInConvergenza && !NonEntrarePiu){
	
	TempoPerArrivareInConvergenza=EntrateInRiallocazione*DimIntervalloRiallocazioniBanda;

	for(i=0;i<NumBuffers_IP; i++){

		TotSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]=0.0;
		TotQuadroSuTuttalaSimDopoConvergenzaPacketLossTramaATM_FlussoIP[i]=0.0;

		MediaPrimaConvergenzaScostamentoPercentualeLoss_IPoATM[i]=TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[i]/EntrateInRiallocazione;
		TotSuTuttalaSimDopoConvergenzaScostamentoPercentualeLoss_IPoATM[i]=0.0;
	}

	//riazzero il conto dei campioni e questo è il 1mo campione
	EntrateInRiallocazione=1;
	
	NonEntrarePiu=true;
}
*/

//END calcoli relativi alle misure medie su tutta la simulazione dell'intervallo di confidenze etc.
////////////////////////////////////////////////
////////////////////////////////////////////////

//mi serve per azzerare i contatori
AzzeraContatoriPerRiprendereIlContoDalRegime();

track TracciaAttuale;
TracciaAttuale.evento = Riallocazione_Ct_ATM;
TracciaAttuale.istante = ClockAttuale + DimIntervalloRiallocazioniBanda;
L.push_back(TracciaAttuale);
L.sort();

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_IdealeLoss_2(double ClockAttuale){

int i;

/////////// "scendo" con la allocazione ideale  ///////////

double BitArrivati_0=Bit_PacchettiIPServiti_0[IndiceletturaPacchettiIPServiti];
double BitArrivati_1=Bit_PacchettiIPServiti_1[IndiceletturaPacchettiIPServiti];

double PlossPrevista_0=PlossoverATM_0[IndiceletturaPacchettiIPServiti];
double PlossPrevista_1=PlossoverATM_1[IndiceletturaPacchettiIPServiti];

Ct_ATM_IPA[0] = (MaxNumConn_IP_AttiveOra-1) * (Init_Bp[0]*( MediaDurataBurst[0]/(MediaDurataBurst[0]+MediaDurataSilenzio[0]) ));
//Ct_ATM_IPA[0]=(1+CellTaxTeorico/100)*Ct_ATM_IPA[0];
Ct_ATM = Ct_ATM_IPA[0]+(BitArrivati_0 / DimIntervalloRiallocazioniBanda) * (PlossPrevista_0-ReferenceIdealeLoss_0);

IndiceletturaPacchettiIPServiti++;
/*diminuendo del 95% il DimIntervalloRiallocazioniBanda succede che esaurisco verso la fine della simulazione i campioni di previsione sul futuro, perciò ricomincio leggendo gli ultimi campioni 
(non posso ricominciare dall'inizio perchè se c'è un cambio di statistiche i campioni alla fine possono essere diversi da quelli misurati all'inizio)
if(IndiceletturaPacchettiIPServiti==NumCampioniPacchettiIPServiti)
	IndiceletturaPacchettiIPServiti=NumCampioniPacchettiIPServiti-10; */

/////////// END "scendo" con la allocazione ideale ///////////

/////////// Stampa situazione ///////////

printf("\n\n\n\t***** RiallocazioneCt_ATM IdealeLoss_2 al clock %lf*****",ClockAttuale);

/*
printf("\n");
for(i=0;i<NumBuffers_IP; i++)
	printf("\nPPerdita_IP[%d]=%e",i,PPerdita_IP[i]);
printf("\n");
printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);
printf("\n");

//Ct_ATM+=ProvaIncremetoCt_ATM;
T_out_ATM =  ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;
printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps", Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps", Ct_IP[1]/((double)(UnMegabit)));

printf("\n\n"); 
//system("pause");
*/

//stampo su file al posto della SpintaGradiente l'e[1]
StampaStato_RiallocazioneCt_ATM(ClockAttuale, -1.0, -1.0, Ct_ATM);

/////////// END Stampa situazione ///////////

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_Ideale(double ClockAttuale){

int i;

/////////// "scendo" con la allocazione ideale  ///////////

//devo conoscere le dimensioni medie x i conti di RCBC successivi
//double DimPacchettoMedio[NumBuffers_IP];
//DimPacchettoMedio[0]=80*UnByte;
//DimPacchettoMedio[1]=4534*UnByte;

double BitArrivati_0=Bit_PacchettiIPServiti_0[IndiceletturaPacchettiIPServiti];
double BitArrivati_1=Bit_PacchettiIPServiti_1[IndiceletturaPacchettiIPServiti];
//BitArrivati_0-=((MaxDimBuffer_ATM-DimBufferAttuale_ATM)/2)*(double)(CellSize)*(double)(DimElementare);
//BitArrivati_1-=((MaxDimBuffer_ATM-DimBufferAttuale_ATM)/2)*(double)(CellSize)*(double)(DimElementare);
double RateIdeale_0=BitArrivati_0 / DimIntervalloRiallocazioniBanda;
double RateIdeale_1=BitArrivati_1 / DimIntervalloRiallocazioniBanda;

//con solo un traffico a livello IP
//Ct_ATM= RateIdeale_0*(1-ReferenceIdealeLoss_0);// + RateIdeale_1;

//con due traffici
double Ct_ATM_0 = RateIdeale_0*(1-ReferenceIdealeLoss_0);
double Ct_ATM_1 = RateIdeale_1*(1-ReferenceIdealeLoss_1);
/*if( Ct_ATM_0 >= Ct_ATM_1 )
	Ct_ATM= Ct_ATM_0;
else
	Ct_ATM= Ct_ATM_1;*/
Ct_ATM= Ct_ATM_0 + Ct_ATM_1;

//un controllo che serve semplicemente ad evitare uno sballo nel calcolo della Ct_ATM in corrispondenza del cambio statistico (anche se non so perchè sbaglia)
if(Ct_ATM<0.8*Ct_ATM__OLD)
	Ct_ATM=Ct_ATM__OLD;
else
	Ct_ATM__OLD=Ct_ATM;

IndiceletturaPacchettiIPServiti++;
/*diminuendo del 95% il DimIntervalloRiallocazioniBanda succede che esaurisco verso la fine della simulazione i campioni di previsione sul futuro, perciò ricomincio leggendo gli ultimi campioni 
(non posso ricominciare dall'inizio perchè se c'è un cambio di statistiche i campioni alla fine possono essere diversi da quelli misurati all'inizio)*/
if(IndiceletturaPacchettiIPServiti==NumCampioniPacchettiIPServiti)
	IndiceletturaPacchettiIPServiti=NumCampioniPacchettiIPServiti-10;

/////////// END "scendo" con la allocazione ideale ///////////

/////////// Stampa situazione ///////////

printf("\n\n\n\t***** RiallocazioneCt_ATM Ideale al clock %lf*****",ClockAttuale);

/*
printf("\n");
for(i=0;i<NumBuffers_IP; i++)
	printf("\nPPerdita_IP[%d]=%e",i,PPerdita_IP[i]);
printf("\n");
printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);
printf("\n");

//Ct_ATM+=ProvaIncremetoCt_ATM;
T_out_ATM =  ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;
printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps", Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps", Ct_IP[1]/((double)(UnMegabit)));

printf("\n\n"); 
//system("pause");
*/

//stampo su file al posto della SpintaGradiente l'e[1]
StampaStato_RiallocazioneCt_ATM(ClockAttuale, -1.0, -1.0, Ct_ATM);

/////////// END Stampa situazione ///////////

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_Bm(double ClockAttuale){

int i;

/////////// "scendo" con la allocazione Bm  ///////////

Ct_ATM_IPA[0] = (MaxNumConn_IP_AttiveOra-1) * (Init_Bp[0]*( MediaDurataBurst[0]/(MediaDurataBurst[0]+MediaDurataSilenzio[0]) ));
//Ct_ATM_IPA[0]=(1+CellTaxTeorico/100)*Ct_ATM_IPA[0];
Ct_ATM = Ct_ATM_IPA[0];

/////////// END "scendo" con la allocazione Bm ///////////

/////////// Stampa situazione ///////////

printf("\n\n\n\t***** RiallocazioneCt_ATM Bm al clock %lf*****",ClockAttuale);

//stampo su file al posto della SpintaGradiente l'e[1]
StampaStato_RiallocazioneCt_ATM(ClockAttuale, -1.0, -1.0, Ct_ATM);

/////////// END Stampa situazione ///////////

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_IPA(double ClockAttuale){

int i;

double SpintaGradiente=0.0;

/////////// StimaDerivata ///////////

if(!Buffer_ATMVuoto){//tengo conto anche del corrente Busy Period

	for(i=0;i<NumBuffers_IP; i++){	
		if(BusyPeriodAttivo[i]){ 
			if(UltimoIstantePerdita_IPoATM[i]!=0.0){

				ContributoIPA_IPoATM[i] += ( UltimoIstantePerdita_IPoATM[i] - InizioBusyPeriod_IPoATM[i] );
				UltimoIstantePerdita_IPoATM[i]=ClockAttuale;
			}
			if(UltimoIstantePerdita_IPoATM_Delay[i]!=0.0){

				ContributoIPA_IPoATM_Delay[i] += ( UltimoIstantePerdita_IPoATM_Delay[i] - InizioBusyPeriod_IPoATM[i] );
				UltimoIstantePerdita_IPoATM_Delay[i]=ClockAttuale;				
			}
		}
		InizioBusyPeriod_IPoATM[i]=ClockAttuale;
	}//end for(i=0;i<NumBuffers_IP; i++)

}// if(!Buffer_ATMVuoto)

// NOTA: il codice di inseguimento della loss che misuro sul buffer IP lo trovi nel backup di questo codice prima della modifica col Delay

// inseguo una LOSS ideale con RCBC
	StimaIPADerivata_IPoATM[0] = (2/DimIntervalloRiallocazioniBanda)*(-ContributoIPA_IPoATM[0])*
										((PacchettiTotaliPersiNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda) 
													- LossRate_IDEALI_0);
	StimaIPADerivata_IPoATM[1] = (2/DimIntervalloRiallocazioniBanda)*(-ContributoIPA_IPoATM[1])*
										((PacchettiTotaliPersiNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda) 
													- LossRate_IDEALI_1);

// inseguo un DELAY ideale con RCBC
	StimaIPADerivata_IPoATM_Delay[0] = (2/DimIntervalloRiallocazioniBanda)*(-ContributoIPA_IPoATM_Delay[0])*
											((PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda) 
													- DelayRate_IDEALI_0);
	StimaIPADerivata_IPoATM_Delay[1] = (2/DimIntervalloRiallocazioniBanda)*(-ContributoIPA_IPoATM_Delay[1])*
											((PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda) 
													- DelayRate_IDEALI_1);

// Nota: l'inseguimento di una loss ideale con PID con solo P, ossia il C2S, lo trovi nel backup di questo codice prima della modifica col Delay

/////////// END StimaDerivata ///////////

/////////// Dimensionamento eta ///////////

//accellero un po la spinta del gradiente
double aiutino = 
0.5*modulo(PacketLossTramaATM_FlussoIP[0]*1e2-ReferenceIdealeLoss_0*1e2);
//0.5*modulo(PacketDelayTramaATM_FlussoIP[0]*0.2*1e2-ReferenceIdealeDelay_0*0.2*1e2);
//pow( ( PacketLossTramaATM_FlussoIP[0]-ReferenceIdealeLoss_0 ), 2);

double AiutinoAlGradientino_0=   aiutino;		
double AiutinoAlGradientino_1=  aiutino;		

StimaIPADerivata_IPoATM[0] = StimaIPADerivata_IPoATM[0] * AiutinoAlGradientino_0;
StimaIPADerivata_IPoATM[1] = StimaIPADerivata_IPoATM[1] * AiutinoAlGradientino_1;

StimaIPADerivata_IPoATM_Delay[0] = StimaIPADerivata_IPoATM_Delay[0] * AiutinoAlGradientino_0;
StimaIPADerivata_IPoATM_Delay[1] = StimaIPADerivata_IPoATM_Delay[1] * AiutinoAlGradientino_1;

/////////// END Dimensionamento eta ///////////

/////////// Scendo Col Gradiente rispetto alla loss solo x la 1ma classe... ///////////

//..., ma attenzione che devo disabilitare l'effetto del video nel buffer L2, dopodichè considero valide le stampe solo della classe 1 (cioè quella del buffer IP 1, tipicamente un VoIP od un on-off)
//Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0];

/////////// END Scendo Col Gradiente rispetto alla loss solo x la 1ma classe ///////////

if(!LavoraSuPDelay){

	/////////// Scendo Col Gradiente rispetto alla loss ///////////

	//Nota: l'esplicitazione degli step di riallocazione ed il loro controllo li trovi nel backup di questo codice prima della modifica col Delay

	if(StimaIPADerivata_IPoATM[1]*StimaIPADerivata_IPoATM[0]>0){   // sono dello stesso segno ...

		if(StimaIPADerivata_IPoATM[1]<0){ // e sono negativi
		
			/*if(-StimaIPADerivata_IPoATM[1]>=-StimaIPADerivata_IPoATM[0]) // spingo in avanti la banda con il modulo + forte
				Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[1];
			else
				Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0];*/
			
			//do la somma
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0] - StimaIPADerivata_IPoATM[1];
		}
		else{ // o sono positivi

			if(-StimaIPADerivata_IPoATM[1]>=-StimaIPADerivata_IPoATM[0]) // spingo indietro la banda con il modulo + debole
				Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0];
			else
				Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[1];
			
			//con la somma
			//Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0] - StimaIPADerivata_IPoATM[1];
		}

	}
	else{ // sono di segno opposto
		
		if(StimaIPADerivata_IPoATM[1]<0) //se e' il video che e' negativo spingo in avanti la banda 
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[1]; //con il video
		else
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM[0];	// altrimenti con il VoIP
	}

	/////////// END Scendo Col Gradiente rispetto alla loss ///////////

}//END if(!LavoraSuPDelay)
else{

/////////// Scendo Col Gradiente rispetto al delay ///////////

//Nota: l'esplicitazione degli step di riallocazione ed il loro controllo li trovi nel backup di questo codice prima della modifica col Delay

if(StimaIPADerivata_IPoATM_Delay[1]*StimaIPADerivata_IPoATM_Delay[0]>0){   // sono dello stesso segno ...

	if(StimaIPADerivata_IPoATM_Delay[1]<0){ // e sono negativi
	
		if(-StimaIPADerivata_IPoATM_Delay[1]>=-StimaIPADerivata_IPoATM_Delay[0]) // spingo in avanti la banda con il modulo + forte
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[1];
		else
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[0];
		
		//do la somma
		//Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[0] - StimaIPADerivata_IPoATM_Delay[1];
	}
	else{ // o sono positivi

		if(-StimaIPADerivata_IPoATM_Delay[1]>=-StimaIPADerivata_IPoATM_Delay[0]) // spingo indietro la banda con il modulo + forte
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[1];
		else
			Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[0];
		
		//do la somma
		//Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[0] - StimaIPADerivata_IPoATM_Delay[1];
	}

}
else{ // sono di segno opposto
	
	if(StimaIPADerivata_IPoATM_Delay[1]<0) //se e' il video che e' negativo spingo in avanti la banda 
		Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[1]; //con il video
	else
		Ct_ATM=Ct_ATM - StimaIPADerivata_IPoATM_Delay[0];	// altrimenti con il VoIP
}
/////////// END Scendo Col Gradiente rispetto al delay ///////////

}//END else dell'if(!LavoraSuPDelay)

printf("\n\n\n\t***** RiallocazioneCt_ATM IPA al clock %lf*****",ClockAttuale);

StampaStato_RiallocazioneCt_ATM(ClockAttuale, ScartoQuadraticoMedio_IPTramaATM, SpintaGradiente, Ct_ATM);

//devo riazzerare i contatori dell'IPA, valutare anche la possibilita' di non azzerarli per la consistenza
for(i=0;i<NumBuffers_IP; i++){
	ContributoIPA_IPoATM[i]=0.0;
	ContributoIPA_IPoATM_Delay[i]=0.0;
}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_BEqBufferless(double ClockAttuale){

int i;

/////////// Rialloco con formula banda equivalente bufferless ///////////

double Media_bitRate_ATM;
double StDev_bitRate_ATM;

/*
RitardoMedio = RitardoTot / CelleServiteTot;

VarianzaRitardo = RitardoTotQuadro + (CelleServiteTot * pow(RitardoMedio,2.0)) - (2.0 * RitardoMedio * RitardoTot);
VarianzaRitardo = VarianzaRitardo / CelleServiteTot;
VarianzaRitardo = sqrt(VarianzaRitardo);
*/

Media_bitRate_ATM=CampioneBitRate_monitoraggioATM / NumCampioni_bitrate;   

StDev_bitRate_ATM = CampioneBitRateQuadro_monitoraggioATM + (NumCampioni_bitrate * pow(Media_bitRate_ATM,2.0)) - (2.0 * Media_bitRate_ATM * CampioneBitRate_monitoraggioATM);
StDev_bitRate_ATM = StDev_bitRate_ATM / NumCampioni_bitrate;
StDev_bitRate_ATM = sqrt(StDev_bitRate_ATM);

double Pi_Greco=3.141592653589793238462643383279502884197169399375;
double a=-2*log(PPerdita_IP_IDEALE)-log(2*Pi_Greco);
a = sqrt(a);

Ct_ATM = Media_bitRate_ATM + a * StDev_bitRate_ATM;
//Ct_ATM	+= Ct_ATM*0.3;

CampioneBitRate_monitoraggioATM = 0.0;   
CampioneBitRateQuadro_monitoraggioATM = 0.0;

/////////// END Rialloco con formula banda equivalente bufferless ///////////

/////////// Stampa situazione ///////////

printf("\n\n\n\t***** RiallocazioneCt_ATM EqB al clock %lf*****",ClockAttuale);

/*for(i=0;i<NumBuffers_IP; i++)
	ScartoQuadraticoMedio_IPTramaATM += pow( ( PacchettiPersiTot[i]-PacchettiTotaliPersiNellaTramaATM_FlussoIP[i] ), 2); 

printf("\n");
for(i=0;i<NumBuffers_IP; i++)
	printf("\nPPerdita_IP[%d]=%e",i,PPerdita_IP[i]);
printf("\n");
printf("\nPacketLossTramaATM=%e",PacketLossTramaATM);
printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);

printf("\nPacketDelayTramaATM_FlussoIP[0]=%e",PacketDelayTramaATM_FlussoIP[0]);
printf("\nPacketDelayTramaATM_FlussoIP[1]=%e",PacketDelayTramaATM_FlussoIP[1]);
printf("\n");

//StampaRisultati();

//Ct_ATM+=ProvaIncremetoCt_ATM;
T_out_ATM =  ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;
printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps", Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps", Ct_IP[1]/((double)(UnMegabit)));
printf("\n\t****** Simulated CellTax: \t%lf ", ( ((Ct_ATM-(Ct_IP[0]+Ct_IP[1]))/(Ct_IP[0]+Ct_IP[1])) * 100) );
printf("\n\n\t****** CellTax teorico: \t%lf ", CellTaxTeorico);
printf("\n\t****** Capacita' del canale ATM secondo CellTax teorico: \t%lf Mbps", Ct_ATMSecondoCellTaxTeorico/((double)(UnMegabit)));
//printf("\n\t****** Tempo per emettere sul canale ATM una cella= %e",T_out_ATM);
*/

//TRUCCO.... metto nel file della SpintaGradiente la Media_bitRate_ATM misurata in real time
StampaStato_RiallocazioneCt_ATM(ClockAttuale, ScartoQuadraticoMedio_IPTramaATM, Media_bitRate_ATM, Ct_ATM);

/////////// END Stampa situazione ///////////

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double RegoleRulex(int N, double bw, double loss, double mean, double var, int B, int MaxB, double Silence, double Bp){

	int out;
	
	//regole di rulex --------- regole con euristica BmBp ---------
	if (mean > 1.722 && loss > 0.020 && Silence <= 7.446) out = 1;
	else if (N > 82 && mean > 0.468 && loss > 0.020 && MaxB > 381) out = 1;
	else if (var <= 0.296 && loss <= 0.000 && bw > 1.751 && B <= 502 && 3.254 < Silence && Silence <= 9.160 && Bp <= 80.225) out = 2;
	else if (loss > 0.020 && bw > 1.751 && Bp <= 71.919) out = 1;
	else if (var > 0.127 && loss > 0.020 && bw <= 2.581 && Silence <= 8.491) out = 1;
	else if (loss <= 0.000 && bw > 4.329 && Silence > 3.691) out = 2;
	else if (loss <= 0.000 && bw > 2.581 && B > 13 && Silence > 6.223 && Bp <= 94.309) out = 2;
	else if (N <= 104 && loss <= 0.002 && bw > 2.581 && 13 < B && B <= 438 && MaxB <= 636 && Silence > 3.691) out = 2;
	else if (var > 0.188 && 0.002 < loss && loss <= 0.020 && bw <= 4.329) out = 0;
	else if (loss > 0.020 && MaxB <= 874 && Silence <= 7.446 && Bp <= 35.620) out = 1;
	else if (var <= 0.158 && loss <= 0.000 && bw > 1.751 && MaxB <= 874 && Silence > 1.865) out = 2;
	else if (N <= 89 && var <= 0.239 && loss <= 0.000 && bw > 1.242 && Silence > 5.557) out = 2;
	else if (N > 89 && mean <= 2.838 && var <= 0.239 && loss <= 0.002 && bw > 1.751 && Silence > 1.865 && Bp <= 66.371) out = 2;
	else if (var > 0.037 && loss > 0.020 && bw <= 1.242 && B > 121) out = 1;
	else if (N > 76 && mean <= 3.625 && var <= 0.373 && loss <= 0.000 && bw > 2.581 && B > 104 && MaxB > 500 && 1.865 < Silence && Silence <= 9.160) out = 2;
	else if (mean <= 0.973 && loss <= 0.000 && bw > 0.410 && B <= 502 && MaxB <= 636 && Silence > 3.691) out = 2;
	else if (N > 82 && mean > 1.233 && loss > 0.002 && bw <= 4.329 && B <= 367 && 95 < MaxB && MaxB <= 765 && Silence <= 4.153) out = 1;
	else if (var > 0.188 && 0.000 < loss && loss <= 0.020 && 43.530 < Bp && Bp <= 71.919) out = 0;
	else if (loss <= 0.002 && bw > 5.375 && B > 43 && Bp <= 97.100) out = 2;
	else if (mean <= 2.343 && var > 0.127 && 0.000 < loss && loss <= 0.020 && bw <= 3.330 && B <= 306 && Silence > 5.557 && Bp > 35.620) out = 0;
	//END regole di rulex

	//applico 3 regole
	if (out==1)
		bw+=(bw*0.15);
	if (out==2)
		bw-=(bw*0.15);
	//END applico 3 regole	

	return bw;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void CalcolaMediaVarianza(){

int i;

double Media_bitRate_ATM;
double StDev_bitRate_ATM;

Media_bitRate_ATM=CampioneBitRate_monitoraggioATM / NumCampioni_bitrate;   

StDev_bitRate_ATM = CampioneBitRateQuadro_monitoraggioATM + (NumCampioni_bitrate * pow(Media_bitRate_ATM,2.0)) - (2.0 * Media_bitRate_ATM * CampioneBitRate_monitoraggioATM);
StDev_bitRate_ATM = StDev_bitRate_ATM / NumCampioni_bitrate;
StDev_bitRate_ATM = sqrt(StDev_bitRate_ATM);

CampioneBitRate_monitoraggioATM = 0.0;   
CampioneBitRateQuadro_monitoraggioATM = 0.0;

Addestramento_Media_bitRate_ATM=Media_bitRate_ATM;
Addestramento_StDev_bitRate_ATM=StDev_bitRate_ATM;

}
////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////// PID ///////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
void RiallocazioneCt_ATM_PID(double ClockAttuale){

int i;

//LO 0 INDICA IL VoIP, L'1 INDICA IL VIDEO

double DeltaU_k[NumBuffers_IP];

/////////// x PID: calcolo loss rate ideale e reale per ogni servizio e l'errore ///////////

// calcolo errore

if(!LavoraSuPDelay){

	e_k[0] = ((PacchettiTotaliPersiNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda) - LossRate_IDEALI_0);
	e_k[1] = ((PacchettiTotaliPersiNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda) - LossRate_IDEALI_1);

}
else{

	e_k[0] = ((PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda) - DelayRate_IDEALI_0);
	e_k[1] = ((PacchettiTotaliInRitardoNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda) - DelayRate_IDEALI_1);
}

/////////// END x PID: calcolo loss rate ideale e reale per ogni servizio e l'errore ///////////

/////////// "Scendo" Col PID ///////////

for(i=0;i<NumBuffers_IP;i++){

	DeltaU_k[i] = (Kp[i]+Ki[i]+Kd[i])*e_k[i] + (-Kp[i]-2*Kd[i])*e_kMeno_1[i] + Kd[i]*e_kMeno_2[i];

	Ct_ATMSeparati[i] = Ct_ATMSeparati[i] + DeltaU_k[i]; 
	//Ct_ATMSeparati[i] = Ct_ATMSeparati[i] + e_k[i]; 
}
/*
printf("\n\n");
printf("\ne_k[0]=%lf", e_k[0]);
printf("\ne_k[1]=%lf", e_k[1]);
printf("\nDeltaU_k[0]=%lf", DeltaU_k[0]);
printf("\nDeltaU_k[1]=%lf", DeltaU_k[1]);
printf("\n");
printf("\nCt_ATMSeparati[0]=%lf Mbps", Ct_ATMSeparati[0]/((double)(UnMegabit)));
printf("\nCt_ATMSeparati[1]=%lf Mbps", Ct_ATMSeparati[1]/((double)(UnMegabit)));
printf("\nLoro somma=%lf Mbps", (Ct_ATMSeparati[0]+Ct_ATMSeparati[1])/((double)(UnMegabit)));
printf("\n\n");
system("pause");
*/

 //se ho due classi di traffico
/*if(Ct_ATMSeparati[0]>Ct_ATMSeparati[1])
	Ct_ATM=Ct_ATMSeparati[0];
else
	Ct_ATM=Ct_ATMSeparati[1];*/
	//in ogni caso, per due classi di traffico bisognerebbe metterci esattamente i confronti fatti x l'RCBC sulle spinte delle due classi (che in questo caso sono il DeltaU_k[i]) 
	//per decidere di quanto aumentare la banda
//oppure do la somma
//Ct_ATM=Ct_ATMSeparati[0]+Ct_ATMSeparati[1];

if(Ct_ATMSeparati[0]>0.0 && Ct_ATMSeparati[1]>0.0){//sono entrambi positivi alzo con la somma come da RCBC papero ToMC
	Ct_ATM=Ct_ATMSeparati[0]+Ct_ATMSeparati[1];
	printf("\nPID do la somma");
	//system("pause");
}
if(Ct_ATMSeparati[0]<=0.0 && Ct_ATMSeparati[1]<=0.0){//sono entrambi negativi abbasso con il piu' debole come da RCBC papero ToMC
	if( modulo(Ct_ATMSeparati[0])>modulo(Ct_ATMSeparati[1]) )
		Ct_ATM=Ct_ATMSeparati[1];
	else
		Ct_ATM=Ct_ATMSeparati[0];
	printf("\nPID abbasso con il piu' debole");
	//system("pause");
}
if( (Ct_ATMSeparati[0]*Ct_ATMSeparati[1])<0.0 ){//uno dei due e' negativo, alzo con quello positivo
	if(Ct_ATMSeparati[0]>0.0)
		Ct_ATM=Ct_ATMSeparati[0];
	if(Ct_ATMSeparati[1]>0.0)
		Ct_ATM=Ct_ATMSeparati[1];
	printf("\nPID uno dei due e' negativo, alzo con quello che devo alzare");
	//system("pause");
}

//se ho una classe di traffico
//Ct_ATM=Ct_ATMSeparati[0];

for(i=0;i<NumBuffers_IP;i++){

	e_kMeno_2[i]=e_kMeno_1[i];
	e_kMeno_1[i]=e_k[i];
}

/*if(e_k[0]>=e_k[1])
	Ct_ATM=Ct_ATM*/

/////////// END "Scendo" Col PID ///////////

/////////// Stampa situazione ///////////

printf("\n\n\n\t***** RiallocazioneCt_ATM PID al clock %lf*****",ClockAttuale);
/*
for(i=0;i<NumBuffers_IP; i++)
	ScartoQuadraticoMedio_IPTramaATM += pow( ( PacchettiPersiTot[i]-PacchettiTotaliPersiNellaTramaATM_FlussoIP[i] ), 2); 

printf("\n");
for(i=0;i<NumBuffers_IP; i++)
	printf("\nPPerdita_IP[%d]=%e",i,PPerdita_IP[i]);
printf("\n");
//printf("\nPacketLossTramaATM=%e",PacketLossTramaATM);
printf("\nPacketLossTramaATM_FlussoIP[0]=%e",PacketLossTramaATM_FlussoIP[0]);
printf("\nPacketLossTramaATM_FlussoIP[1]=%e",PacketLossTramaATM_FlussoIP[1]);
printf("\ne[0]=%lf Kbps", e_k[0]/((double)(UnKilobit)));
printf("\ne[1]=%lf Kbps", e_k[1]/((double)(UnKilobit)));

printf("\nLossRateNellaTramaATM_FlussoIP[0]=%lf Kbps", (PacchettiTotaliPersiNellaTramaATM_FlussoIP[0]*DimPacchettoMedio[0]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nLossRateNellaTramaATM_FlussoIP[1]=%lf Kbps", (PacchettiTotaliPersiNellaTramaATM_FlussoIP[1]*DimPacchettoMedio[1]/DimIntervalloRiallocazioniBanda)/((double)(UnKilobit)));
printf("\nLossRate_IDEALE_0=%lf Kbps", LossRate_IDEALI_0/((double)(UnKilobit)));
printf("\nLossRate_IDEALE_1=%lf Kbps", LossRate_IDEALI_1/((double)(UnKilobit)));
printf("\n");

//StampaRisultati();

//Ct_ATM+=ProvaIncremetoCt_ATM;
T_out_ATM =  ((double)(CellSize)*(double)(DimElementare)) / Ct_ATM;
printf("\n\t****** Capacita' del canale ATM: \t%lf Mbps", Ct_ATM/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 0: \t%lf Mbps", Ct_IP[0]/((double)(UnMegabit)));
printf("\n\t****** Capacita' del canale IP 1: \t%lf Mbps", Ct_IP[1]/((double)(UnMegabit)));
//printf("\n\t****** Simulated CellTax: \t%lf ", ( ((Ct_ATM-(Ct_IP[0]+Ct_IP[1]))/(Ct_IP[0]+Ct_IP[1])) * 100) );
//printf("\n\n\t****** CellTax teorico: \t%lf ", CellTaxTeorico);
//printf("\n\t****** Capacita' del canale ATM secondo CellTax teorico: \t%lf Mbps", Ct_ATMSecondoCellTaxTeorico/((double)(UnMegabit)));
//printf("\n\t****** Tempo per emettere sul canale ATM una cella= %e",T_out_ATM);

printf("\n\n"); 
//system("pause");
*/

//stampo su file al posto della SpintaGradiente l'e[1]
StampaStato_RiallocazioneCt_ATM(ClockAttuale, ScartoQuadraticoMedio_IPTramaATM, e_k[1], Ct_ATM);

/////////// END Stampa situazione ///////////

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void StampaStato_RiallocazioneCt_ATM(double ClockAttuale, double ScartoQuadraticoMedio_IPTramaATM, double SpintaGradiente, double Ct_ATM){

	/*

	if(PlossNellaTramaATM_0>0){//controllo x eliminare stampe sballate dovute alla riallocazione proprio sovrapposta al campio statistiche

		fprintf(Time,"%lf\n",ClockAttuale);
		fprintf(PPerditaFlussoIP_0,"%e\n",PPerdita_IP[0]);
		fprintf(PPerditaFlussoIP_1,"%e\n",PPerdita_IP[1]);
		fprintf(ScartoQuadraticoMedio,"%e\n",ScartoQuadraticoMedio_IPTramaATM);
		fprintf(ValoreSpintaGradiente,"%e\n",SpintaGradiente);
		fprintf(ValoreCt_ATM,"%lf\n",Ct_ATM / ((double)UnMegabit));

		fprintf(PlossNellaTramaATM,"%e\n",PacketLossTramaATM);
		
		if(!LavoraSuPDelay){//stampo su file la metrica probabilistica di interesse Loss o Delay anche se tengo sempre il nome relativo alla loss
			fprintf(PlossNellaTramaATM_0,"%e\n",PacketLossTramaATM_FlussoIP[0]);
			fprintf(PlossNellaTramaATM_1,"%e\n",PacketLossTramaATM_FlussoIP[1]);
		}
		else{
			fprintf(PlossNellaTramaATM_0,"%e\n",PacketDelayTramaATM_FlussoIP[0]);
			fprintf(PlossNellaTramaATM_1,"%e\n",PacketDelayTramaATM_FlussoIP[1]);
		}
				
		fprintf(SimulatedCellTax,"%e\n", ( ((Ct_ATM-(Ct_IP[0]+Ct_IP[1]))/(Ct_IP[0]+Ct_IP[1])) * 100));
		fprintf(ValoreCt_IP_0,"%lf\n", Ct_IP[0] / ((double)UnMegabit));
		fprintf(ValoreCt_IP_1,"%lf\n", Ct_IP[1] / ((double)UnMegabit));

		if(TecnicaRiallocazione==0 || TecnicaRiallocazione==4){
			
			//if(Bit_PacchettiServitiTot[0]==0.0 || Bit_PacchettiServitiTot[1]==0.0)
			if(Bit_PacchettiServitiTot[0]==0.0 ){
				printf("\n\nBit_PacchettiServitiTot di 0 od 1 a zero!");
				//system("pause");
			}
			
			fprintf(Bit_PackIPServiti_0,"%e\n",Bit_PacchettiServitiTot[0]);
			fprintf(Bit_PackIPServiti_1,"%e\n",Bit_PacchettiServitiTot[1]);
		}//END if(TecnicaRiallocazione!=3)
	
	} //END controllo x eliminare stampe sballate dovute alla riallocazione proprio sovrapposta al campio statistiche

	*/

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void CambioStatistiche(double ClockAttuale){

int i;
track TracciaAttuale;

printf("\n\n\nSono in CambioStatistiche al tempo %lf",ClockAttuale);
printf("\nCt_IP prima=%lf",Ct_IP[0]/(double)UnMegabit);
printf("\nCt_ATM prima=%lf",Ct_ATM/(double)UnMegabit);

if(GeneraConBurst){

int MaxNumConn_IP_AttivePrima=MaxNumConn_IP_AttiveOra;

if(ClockAttuale==(1*DurataSimulazione/5)){
	
	MaxNumConn_IP_AttiveOra=MaxNumConn_IP-40;

	/*Ct_IP[0]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[0]*0.452);
	Ct_IP[1]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[1]*0.452);*/
	
	for(i=MaxNumConn_IP_AttivePrima; i<MaxNumConn_IP_AttiveOra; i++){
        TracciaAttuale.evento = Inizio_Burst;
        TracciaAttuale.istante = ClockAttuale;
        //TracciaAttuale.istante = ClockAttuale + GenDurataSilenzio(i);
        //TracciaAttuale.istante = QuandoPrimoBurst[i];
        TracciaAttuale.Connessione = i;
        printf("\nPer la connessione %d primo InizioBurst programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
		//printf("\nClock=%lf next burst per la conn %d al t=%lf",0.0,TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
	}
}
if(ClockAttuale==(2*DurataSimulazione/5)){

	MaxNumConn_IP_AttiveOra=MaxNumConn_IP-30;
	
	/*Ct_IP[0]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[0]*0.442);
	Ct_IP[1]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[1]*0.442);*/

	for(i=MaxNumConn_IP_AttivePrima; i<MaxNumConn_IP_AttiveOra; i++){
        TracciaAttuale.evento = Inizio_Burst;
        TracciaAttuale.istante = ClockAttuale;
        //TracciaAttuale.istante = ClockAttuale + GenDurataSilenzio(i);
        //TracciaAttuale.istante = QuandoPrimoBurst[i];
        TracciaAttuale.Connessione = i;
        printf("\nPer la connessione %d primo InizioBurst programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
		//printf("\nClock=%lf next burst per la conn %d al t=%lf",0.0,TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
	}
}
if(ClockAttuale==(3*DurataSimulazione/5)){
	
	MaxNumConn_IP_AttiveOra=MaxNumConn_IP-20;
	
	/*Ct_IP[0]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[0]*0.442);
	Ct_IP[1]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[1]*0.442);*/

	for(i=MaxNumConn_IP_AttivePrima; i<MaxNumConn_IP_AttiveOra; i++){
        TracciaAttuale.evento = Inizio_Burst;
        TracciaAttuale.istante = ClockAttuale;
        //TracciaAttuale.istante = ClockAttuale + GenDurataSilenzio(i);
        //TracciaAttuale.istante = QuandoPrimoBurst[i];
        TracciaAttuale.Connessione = i;
        printf("\nPer la connessione %d primo InizioBurst programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
		//printf("\nClock=%lf next burst per la conn %d al t=%lf",0.0,TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
	}
}
if(ClockAttuale==(4*DurataSimulazione/5)){
	
	MaxNumConn_IP_AttiveOra=MaxNumConn_IP-10;
	
	/*Ct_IP[0]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[0]*0.442);
	Ct_IP[1]=(MaxNumConn_IP_AttiveOra/2)*(Init_Bp[1]*0.442);*/

	for(i=MaxNumConn_IP_AttivePrima; i<MaxNumConn_IP_AttiveOra; i++){
        TracciaAttuale.evento = Inizio_Burst;
        TracciaAttuale.istante = ClockAttuale;
        //TracciaAttuale.istante = ClockAttuale + GenDurataSilenzio(i);
        //TracciaAttuale.istante = QuandoPrimoBurst[i];
        TracciaAttuale.Connessione = i;
        printf("\nPer la connessione %d primo InizioBurst programmato al  tempo:\t%lf",TracciaAttuale.Connessione,TracciaAttuale.istante);
		//printf("\nClock=%lf next burst per la conn %d al t=%lf",0.0,TracciaAttuale.Connessione,TracciaAttuale.istante);
        L.push_back(TracciaAttuale);
	}
}

////////////////////// 20 luglio 2011, mi accorgo ora che al cambio di statistiche la Bp era da aggiornare perche' la funzione CheBufferIPCorrente(i) lavora solo fino al MaxNumConn_IP_AttiveOra
//xcio' nell'inizializza iniziale sebbe il ciclo for fosse fino a MaxNumConn_IP metteva giusto la Bp solo alle MaxNumConn_IP_AttiveOra iniziali cioe' 70
for(i=0; i<MaxNumConn_IP; i++)
		Bp[i]=Init_Bp[CheBufferIPCorrente(i)];

for(i=0; i<MaxNumConn_IP; i++)
	T_in_IP[i]=( (double)(DimPacchetto) ) / Bp[i];

/*printf("\n\n");
for(i=0; i<MaxNumConn_IP_AttiveOra; i++){
	printf("\nBp[%d]: %lf",i,Bp[i]);
	printf("\nT_in_IP[%d]: %lf",i,T_in_IP[i]);
}*/
////////////////////// END 20 luglio 2011

//T_in_IP=( (double)(DimPacchetto) ) / Bp;
for(i=0;i<NumBuffers_IP; i++)
	T_out_IP[i]=( (double)(DimPacchetto) ) / Ct_IP[i];

InizializzaT_in_ATM();

/*Ct_ATM_IPA[0] = (MaxNumConn_IP_AttiveOra-1) * (Init_Bp[0]*( MediaDurataBurst[0]/(MediaDurataBurst[0]+MediaDurataSilenzio[0]) ));
Ct_ATM_IPA[0]=(1+CellTaxTeorico/100)*Ct_ATM_IPA[0];
Ct_ATM_IPA[1] = 200.0* UnKilobit;

Ct_ATM = Ct_ATM_IPA[0] + Ct_ATM_IPA[1] ;*/

T_out_ATM = (((double)(CellSize))*(double)(DimElementare)) / Ct_ATM;   //tempo per emettere sul canale ATM una cella ATM

}// END del if(GeneraConBurst)
/*else{

	Ct_IP[0]=19.0*UnMegabit;
	CtInSorgente_IP=2.0*UnMegabit;
	T_out_IP[0] = ( (double)(DimElementare) ) / Ct_IP[0]; //tempo per emettere sul canale IP una cella elementare
	T_out_IP[1] = ( (double)(DimElementare) ) / Ct_IP[1]; //tempo per emettere sul canale IP una cella elementare
//	T_in_IP = ( (double)(DimElementare) ) / (CtInSorgente_IP); //tempo che la sorgente impiega per inserire sul canale una cella elementare

	Ct_ATM=15.0*UnMegabit;
	T_out_ATM =  (((double)(CellSize))*(double)(DimElementare)) / Ct_ATM;   //tempo per emettere sul canale ATM una cella ATM
//	T_in_ATM =  T_in_IP*((double)(CellSize));

}// END del else (GeneraConBurst)*/

/*printf("\n\nCt_IP[0] ora=%lf",Ct_IP[0]/(double)UnMegabit);
printf("\n\nCt_IP[1] ora=%lf",Ct_IP[1]/(double)UnMegabit);
printf("\nCt_ATM ora=%lf",Ct_ATM/(double)UnMegabit);*/

AggiornaContatoriTot_IP();
AggiornaContatoriTot_ATM();
AzzeraContatoriPerRiprendereIlContoDalRegime();

L.sort();//system("pause");

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void Monitorbitrate(double ClockAttuale){

track TracciaAttuale;

double CampioneBitRateAttuale=BitArrivati_monitoraggioATM / DimIntervalloMonitor_bitrate;

CampioneBitRate_monitoraggioATM += CampioneBitRateAttuale;   
CampioneBitRateQuadro_monitoraggioATM += pow(CampioneBitRateAttuale, 2.0);

BitArrivati_monitoraggioATM=0.0;

TracciaAttuale.evento = Monitor_bitrate;
TracciaAttuale.istante = ClockAttuale+DimIntervalloMonitor_bitrate;
L.push_back(TracciaAttuale);

L.sort();
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InizioBurst(int NumConn, double ClockAttuale){

DurataAttualeBurst[NumConn] = GenDurataBurst(NumConn);
IstanteFineAttualeBurst[NumConn] = ClockAttuale + DurataAttualeBurst[NumConn];

//printf("\nClock=%lf programmato burst per la conn %d con termine al t=%lf", ClockAttuale, NumConn, IstanteFineAttualeBurst[NumConn] );

NumConnAttive++;
ArrivoPacchetto(NumConn,ClockAttuale);
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double GenDurataBurst(int NumConn){

double n,EX;

int FlussoIPCorrente=CheBufferIPCorrente(NumConn);

	do
	{
		n = NumeroRand(1.0);
	}
	while (n==0.0);

        if(GeneraSecondoClaudio)
                EX = deltaMediaDurataBurst[FlussoIPCorrente] * (double)pow(n, 1/(-alfa));

        if(GeneraSecondoMario)
				EX = MediaDurataBurst[FlussoIPCorrente] * (-log(n));

//EX = ( (int)(EX/T_in_IP) ) * T_in_IP;

//printf("\n%lf",EX);

//NumBurstGenerati++;
//DurataEffettivaDeiBurst += EX;

return EX;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
double GenDurataSilenzio(int NumConn){

double a,EX;

int FlussoIPCorrente=CheBufferIPCorrente(NumConn);

    do
	{
		a = NumeroRand(1.0);
	}
	while (a==0.0);

        if(GeneraSecondoClaudio)
				EX = deltaMediaDurataSilenzio[FlussoIPCorrente] * (double)pow(a, 1/(-alfa));

        if(GeneraSecondoMario)
				EX = MediaDurataSilenzio[FlussoIPCorrente] * -(log(a));
				
//EX = ( (int)(EX/T_in_IP) ) * T_in_IP;

//NumSilenziGenerati++;
//DurataEffettivaDeiSilenzi += EX;

//printf("\n%lf",EX);

return EX;
}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InizializzaTracciatoVideo(){ 
int i;


for(i=0;i<NumCampioniTracciatoVideo;++i){
	fscanf(FileTracciatoVideo_Pacchetti,"%lf\n", &CampioniTracciatoVideo_Pacchetti[i]);
	//printf("\nCampioniTracciatoVideo_Pacchetti[%d]=%lf",i,CampioniTracciatoVideo_Pacchetti[i]);
}

for(i=0;i<NumCampioniTracciatoVideo;++i){
	fscanf(FileTracciatoVideo_Tempi,"%lf\n", &CampioniTracciatoVideo_Tempi[i]);
	//printf("\nCampioniTracciatoVideo_Tempi[%d]=%lf",i,CampioniTracciatoVideo_Tempi[i]);
}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InizializzaTracciatoPacchettiIPServiti(){ 

	int i;

	for(i=0;i<NumCampioniPacchettiIPServiti;++i){

		fscanf(Bit_PackIPServiti_0,"%lf\n", &Bit_PacchettiIPServiti_0[i]);
		fscanf(Bit_PackIPServiti_1,"%lf\n", &Bit_PacchettiIPServiti_1[i]);

	}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
void InizializzaTracciatoPloss(){ 

	int i;

	for(i=0;i<NumCampioniPacchettiIPServiti;++i){

		fscanf(PlossNellaTramaATM_0,"%lf\n", &PlossoverATM_0[i]);
		fscanf(PlossNellaTramaATM_1,"%lf\n", &PlossoverATM_1[i]);

	}

}
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
int ConvertiAdInteroSuperiore(double num){

int out;

if( (num-(int)num) > 0.0 )
	out=(int)num+1;
else
	out=(int) num;

return out;
}
////////////////////////////////////////////////////////////////////

double modulo(double num){

	if(num<0)
		return (-num);
	else
		return num;
}


////////////////////////////////////////////////////////////////////
void SkipSobol(FILE** FileSobol_ref, int skip) {
  int i, j;
  double fake;

  for (j = 0; j < skip; j++) {
    for (i = 0; i < 6; i++) {
      if (i < 5)
        fscanf(*FileSobol_ref, "%lf\t", &fake);
      else
        fscanf(*FileSobol_ref, "%lf\n", &fake);
    }
  }
}




