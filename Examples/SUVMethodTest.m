function [] = SUVMethodTest()
%SUVMETHODTEST
paramNames = {'methodTest','dim','paramNote','MS_maxIter','includePhase','searchAlgorithm','DESCENT_Gamma','nIC','MS_nLM'}
paramValues = {true,2,'methodTest',30,false,"Stochastic Unit Vectors",1,4,4}
data = Main(paramNames,paramValues);

end

