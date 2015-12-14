
package ChemProcSim

// HACK All streams are upper bounded which is ugly. It would be nice to 
// bound streams on a composite unit level
  type FlowRate = Real(min=0, max=10.0, start=0.1, quantity = "molar flowrate", unit = "mol/s");
//  type Pressure = Real(min = 0, start=1, quantity = "pressure", unit = "Pa");
  type EnthalpyFlow = Real(quantity = "enthalpy flowrate", unit = "J/s");
  type MolarEnthalpy = Real(quantity = "molar enthalpy", unit = "J/mol");
  type MolFraction = Real(min=0, max=1, unit="mol/mol");

  constant Integer C=2;  // FIXME


  connector Stream
    FlowRate f[C];
//    Pressure p;
    EnthalpyFlow H;
  end Stream;

  
  connector Inlet  = input  Stream;
  connector Outlet = output Stream;


  class Source
    Outlet outlet;
  end Source;


  class Sink
    Inlet inlet;
  end Sink;


//   partial class PressureChange
//     Pressure difference;
//   end PressureChange;
// 
// 
//   class NoPressureChange
//     extends PressureChange;
//     equation
//       difference = 0.0;
//   end NoPressureChange;


  partial class HeatExchanged
    EnthalpyFlow heat;
  end HeatExchanged;


  class NoHeatExchanged
    extends HeatExchanged;
    equation
      heat = 0.0;
  end NoHeatExchanged;


  partial class Reaction
    FlowRate     rate[C];
    EnthalpyFlow heat;
  end Reaction;


  class NoReaction
    extends Reaction;
    equation
      rate = fill(0.0, C);
      heat = 0.0;
  end NoReaction;


  partial class UnitOp

    replaceable class RateEquations    = NoReaction       constrainedby Reaction;
    replaceable class ExchangedHeat    = NoHeatExchanged  constrainedby HeatExchanged;
//    replaceable class ChangeInPressure = NoPressureChange constrainedby PressureChange;


    RateEquations    reaction ;
    ExchangedHeat    exchanged;
//    ChangeInPressure pressure ;

    parameter Integer nInlet;
    parameter Integer nOutlet;
    
    Inlet[nInlet] inlet;
    Outlet[nOutlet] outlet;

//    Pressure p;
    
    equation
      
        for i in 1:C loop
            sum( inlet[j].f[i] for j in 1:nInlet) =
            sum(outlet[k].f[i] for k in 1:nOutlet) + reaction.rate[i];
        end for;

        sum(inlet[j].H for j in 1:nInlet) = 
        sum(outlet[k].H for k in 1:nOutlet) + reaction.heat + exchanged.heat;

//        p = min(inlet[j].p for j in 1:nInlet) + pressure.difference;
    
//        for k in 1:nOutlet loop
//            outlet[k].p = p;
//        end for;

  end UnitOp;


  partial class SISO  // TODO Why cannot this be inlined?
    extends UnitOp(nInlet=1, nOutlet=1);
  end SISO;


  class HeatExchanger
    extends SISO;
    redeclare class ExchangedHeat = HeatExchanged;
  end HeatExchanger;


//   class PressureChanger
//     extends SISO;
//     redeclare class ChangeInPressure = PressureChange;
//   end PressureChanger;


  partial class SIDO
    extends UnitOp(nInlet=1, nOutlet=2);
  end SIDO;


  class Divider
    extends SIDO;
    Real zeta(min=0, max=1);
    equation
      outlet[1].f = zeta*inlet[1].f;
      outlet[1].H = zeta*inlet[1].H;
  end Divider;
  

  partial class FlashBase
    extends SIDO;
    FlowRate L(min=0.1, max=10.0), V(min=0.1, max=10.0); // TODO These should be non-negative too
    // TODO Introduce classes: Phase equilibrium, Enthalpy model
    MolFraction[C] x, y;
    MolarEnthalpy hV, hL; 
    //Real T; // TODO Not sure the best place is here
              //  if not needed (e.g. ideal binary mixtures) then messes up the DoF
    equation
      V = sum(outlet[1].f);
      y = outlet[1].f / V;
      L = sum(outlet[2].f);
      x = outlet[2].f / L;
      outlet[1].H = V*hV;
      outlet[2].H = L*hL;
  end FlashBase;
  
  
  partial class NoVaporFlashBase
    extends SISO;
    FlowRate L(min=0.1, max=10.0);
    MolFraction[C] x;
    MolarEnthalpy hL;
    equation
      L = sum(outlet[1].f);
      x = outlet[1].f / L;
      outlet[1].H = L*hL;
  end FlashBase;


  class ReactiveFlash
    extends FlashBase;
    redeclare class RateEquations = Reaction;
  end ReactiveFlash;


  class Mixer
    extends UnitOp(nInlet=2, nOutlet=1);
  end Mixer;

  class Mixer3
    extends UnitOp(nInlet=3, nOutlet=1);
  end Mixer3;

  // Composite units
  //====================================================================

  class DummyFeed // HACK
    Outlet outlet(f = fill(0.0, C), H = 0);  // FIXME Removed pressure
  end DummyFeed;

  class VLEStage
    
    Inlet[3] inlet;   // L, F, V
    Outlet[2] outlet; // V, L

    replaceable class FlashUnit = FlashBase;
    FlashUnit flash;

    Mixer3 mixer;

    equation
      connect(inlet, mixer.inlet);
      connect(mixer.outlet[1], flash.inlet[1]);
      connect(flash.outlet, outlet);

  end VLEStage;


  class VLEcascade

    Inlet[3]  inlet; // L, F, V
    Outlet[2] outlet; // V, L

    DummyFeed dummyFeed; // FIXME violates the 1 connection / connector rule...
    
    replaceable class Flash = FlashBase;

    class StageClass = VLEStage(redeclare class FlashUnit = Flash);
    
    parameter Integer nStages;
    parameter Integer feedStage;

    StageClass[nStages] stages;

    equation

      connect(stages[1].inlet[1],  inlet[1]);
      connect(stages[1].outlet[1], outlet[1]);
      connect(stages[1].inlet[3], stages[2].outlet[1]);

      for i in 2:nStages-1 loop  // FIXME Magic numbers, enum?
        connect(stages[i].inlet[1], stages[i-1].outlet[2]); // L
        connect(stages[i].inlet[3], stages[i+1].outlet[1]); // V
      end for;

      // FIXME Code triplication
      for i in 1:feedStage-1 loop
        connect(stages[i].inlet[2], dummyFeed.outlet);
      end for;

      connect(stages[feedStage].inlet[2], inlet[2]);

      for i in feedStage+1:nStages loop
        connect(stages[i].inlet[2], dummyFeed.outlet);
      end for;

      connect(stages[nStages].inlet[1], stages[nStages-1].outlet[2]);
      connect(stages[nStages].inlet[3],  inlet[3]);
      connect(stages[nStages].outlet[2], outlet[2]);

  end VLEcascade;


  class TotalCondenser
    Inlet inlet;
    Outlet outlet[2]; // 1: distillate, 2: reflux
    replaceable class FlashUnit = NoVaporFlashBase;
    FlashUnit     noVaporFlash;
    HeatExchanger heatExchanger;
    Divider       divider;
    equation
      connect(inlet, heatExchanger.inlet[1]);
      connect(heatExchanger.outlet[1], noVaporFlash.inlet[1]);
      connect(noVaporFlash.outlet[1], divider.inlet[1]);
      connect(divider.outlet[1], outlet[1]);
      connect(divider.outlet[2], outlet[2]);
  end TotalCondenser;

  // TODO How can I connect a feed to a cascade on the GUI

end ChemProcSim;

model JacobsenTest

    import ChemProcSim.*;

    parameter Real[C] MolWeight = { 32.04, 60.10};

    class IdealBinary
      extends FlashBase;
      parameter Real alpha = 3.55;
      equation
        y[1] = alpha*x[1]/(1.0+(alpha-1.0)*x[1]);
        hV = 0.1349*exp(-3.98*x[1])+0.4397*exp(-0.088*x[1]);
        hL = 0.1667*exp(-1.087*x[1]);
    end IdealBinary;

    class Reboiler
      extends IdealBinary;
      redeclare class ExchangedHeat = HeatExchanged;
    end Reboiler;

    class BubblePointFlash
      extends NoVaporFlashBase;
      equation
        hL = 0.1667*exp(-1.087*x[1]);    
    end BubblePointFlash;
    
    class Condenser
      extends TotalCondenser(redeclare class FlashUnit = BubblePointFlash);
    end Condenser;

    //==================================================================

    Source        feed(outlet.f = {0.5, 0.5}, outlet.H = 0.1667*exp(-1.087*0.5) );  // FIXME Removed pressure
    Condenser     condenser;
    VLEcascade    cascade(redeclare class Flash = IdealBinary, nStages=50, feedStage=25);
    Reboiler      reboiler;
    Sink          distillateSink;
    Sink          liquidSink;

    equation

      connect(cascade.outlet[1], condenser.inlet);
      connect(condenser.outlet[1], distillateSink.inlet);
      connect(condenser.outlet[2], cascade.inlet[1]);

      sum(condenser.outlet[2].f[i]*MolWeight[i] for i in 1:C) = 96.0;

      connect(feed.outlet, cascade.inlet[2]);

      reboiler.V = 3.0;

      connect(cascade.inlet[3], reboiler.outlet[1]);
      connect(cascade.outlet[2],reboiler.inlet[1]);
      connect(reboiler.outlet[2], liquidSink.inlet);

end JacobsenTest;

