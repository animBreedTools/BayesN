###compatible with 1.4.1

using DataFrames
using Distributions
using DelimitedFiles
using LinearAlgebra
using CSV
using Printf

function bayesPR_selReg(genoTrain, phenoTrain, snpInfo, chrs,locusID, fixedRegSize, priorPi, estPi, varGenotypic, varResidual, chainLength, burnIn, outputFreq, onScreen)
    SNPgroups   = prepRegionData(snpInfo, chrs, locusID, fixedRegSize)
    these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations
    nRegions    = length(SNPgroups)
    println("number of regions: ", nRegions)
    dfEffectVar = 4.0
    dfRes       = 4.0
    X           = convert(Array{Float64}, genoTrain)
    println("X is this size", size(X))
    y           = convert(Array{Float64}, phenoTrain)
    println("y is this size", size(y))
    nTraits, nRecords , nMarkers   = size(y,2), size(y,1), size(X,2)
    fileControlSt(fixedRegSize)
    p           = mean(X,dims=1)./2.0
    sum2pq      = sum(2*(1 .- p).*p)
    if varGenotypic==0.0
        varBeta      = fill(0.0005, nRegions)
        else varBeta = fill(varGenotypic/sum2pq, nRegions)
    end
    if varResidual==0.0
        varResidual  = 0.0005
    end
    scaleVar        = varBeta[1]*(dfEffectVar-2.0)/dfEffectVar
    νS_β            = scaleVar*dfEffectVar
    df_β            = dfEffectVar
    scaleRes        = varResidual*(dfRes-2.0)/dfRes
    νS_e            = scaleRes*dfRes
    df_e            = dfRes
    logPiD          = log(priorPi)
    logPiDComp      = log(1-priorPi)
    tempBetaVec     = zeros(Float64,nMarkers) #initial values as "0"
    μ               = mean(y)
    X              .-= ones(Float64,nRecords)*2p
    xpx             = diag(X'X)
    ycorr           = y .- μ
    #MCMC starts here
    for iter in 1:chainLength
        #sample residual variance
        varE = sampleVarE(νS_e,ycorr,df_e,nRecords)
        #sample intercept
        ycorr    .+= μ
        rhs      = sum(ycorr)
        invLhs   = 1.0/nRecords
        meanMu   = rhs*invLhs
        μ        = rand(Normal(meanMu,sqrt(invLhs*varE)))
        ycorr    .-= μ
        regCounter = 0
        for r in 1:nRegions
            theseLoci = SNPgroups[r]
            regionSize = length(theseLoci)
            λ_r = varE/varBeta[r]
            ##select region
            ycorr .+= view(X,:,theseLoci)*tempBetaVec[theseLoci] 
            rhsReg = view(X,:,theseLoci)'*ycorr
            regionXpX = Matrix(Diagonal(xpx[theseLoci])) #can be done initially
            varD0 = regionXpX.*varE
            varD1 = ^(regionXpX,2).*varBeta[r] + varD0
            logD0 = -(0.5)*(regionSize*log(det(varD0)) + rhsReg'*inv(varD0)*rhsReg) + logPiD
            logD1 = -(0.5)*(regionSize*log(det(varD1)) + rhsReg'*inv(varD1)*rhsReg) + logPiDComp
            probD1 = 1.0/(1.0 + exp(logD0-logD1))
            println(probD1)
            ycorr .-= view(X,:,theseLoci)*tempBetaVec[theseLoci]
            if probD1 < rand()
                println("region $r fitted")
                regCounter += 1
                for l in theseLoci::UnitRange{Int64}
                   BLAS.axpy!(tempBetaVec[l], view(X,:,l), ycorr)
                   rhs = view(X,:,l)'*ycorr
                   lhs = xpx[l] + λ_r
                   meanBeta = lhs\rhs
                   tempBetaVec[l] = sampleBeta(meanBeta, lhs, varE)
                   BLAS.axpy!(-1*tempBetaVec[l], view(X,:,l), ycorr)
                end
            else
                println("region $r NOT fitted")
                tempBetaVec[theseLoci] .= 0
            end
            varBeta[r] = sampleVarBeta(νS_β,tempBetaVec[theseLoci],df_β,regionSize)
        end
        println("fitted regions: $regCounter, prop fitted regions: $(regCounter/nRegions)")
        outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    end
#    betaFromFile =  readcsv(pwd()"/betaOut",header=false)
#    print("read $(size(betaFromFile,1)) samples for $(size(betaFromFile,2)) markers from betaOut \n")
#    bayesRegOut = mean(betaFromFile,1)
#    @printf("acc %.6f \n", cor(y,X*bayesRegOut')[1])
end

###now genoTrain is excluded, and it takes locusIDs, and map to trim
function prepRegionData(userMapData,chrs,locusID,fixedRegSize)
    accRegion = 0
    accRegionVec = [0]
    SNPgroups = []
    headMap = [:row, :snpID, :snpOrder ,:chrID, :pos]
    #for Ana's map
    mapData = userMapData
    #
    rename!(mapData, headMap)     ############## names(mapData) removed ############
    print(mapData[1:5,:])
    print(mapData[1:10,:])
    ###
    mapData = mapData[mapData[!,:chrID] .<= chrs,:]
    # if first col in genoTrain is ID
    # I find cols that are in mapData (<chrs), and select those
#    usedLoci = intersect(Symbol.(locusID),Symbol.(mapData[!,:snpID]))
#    mapData = mapData[[findall(usedLoci[i].==Symbol.(mapData[!,:snpID]))[] for i in 1:length(usedLoci)],:] #trim map data  ######## find -> findall #########
#    println([findall(usedLoci[i].==Symbol.(mapData[!,:snpID]))[] for i in 1:length(usedLoci)])    ######## just added to check #######  
#    totLoci = length(usedLoci) # first col is ID

    tempLocusID = DataFrame()
    tempLocusID.snpID = locusID
    mapData = rightjoin(mapData,tempLocusID,on=:snpID)
    totLoci = size(mapData,1)
    
    println("totalLoci in MAP: $totLoci")
    snpInfoFinal = DataFrame([Vector{Any}(undef, 0) for i = 1:3])
    if fixedRegSize==99
        println("fixedRedSize $fixedRegSize")
        snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
        accRegion    = length(unique(mapData[!,:chrID]))
        elseif fixedRegSize==9999
            snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
            snpInfoFinal[!,:chrID]  .= 1 #was "1"
            accRegion    = 1
        else
        for c in 1:chrs
            thisChr = mapData[mapData[!,:chrID] .== c,:]
            totLociChr = size(thisChr,1)
            TotRegions = ceil(Int,totLociChr/fixedRegSize)
            accRegion += TotRegions
            push!(accRegionVec, accRegion)
            tempGroups = sort(repeat(collect(accRegionVec[c]+1:accRegionVec[c+1]),fixedRegSize))
            snpInfo = DataFrame(Any, length(tempGroups), 3)
            snpInfo[1:totLociChr,1] = collect(1:totLociChr)
            snpInfo[1:totLociChr,2] = thisChr[!,:snpID]
            snpInfo[:,3] = tempGroups
            snpInfo = snpInfo[[isassigned(snpInfo[:,1],i) for i in 1:size(snpInfo,1)],:]
            snpInfoFinal = vcat(snpInfoFinal,snpInfo)
            @printf("chr %.0f has %.0f groups \n", c, TotRegions)
            println(by(snpInfo, :x3, nrow)[:,2])
        end
        end  #ends if control flow
#    print(snpInfoFinal)
    CSV.write("snpInfo",snpInfoFinal,header=false)
    for g in 1:accRegion
        push!(SNPgroups,searchsorted(snpInfoFinal[:,3], g))
    end
    return SNPgroups #, genoX
end

function outputControlSt(onScreen,iter,these2Keep,X,tempBetaVec,μ,varBeta,varE,fixedRegSize)
    if iter in these2Keep
        out0 = open(pwd()*"/muOut$fixedRegSize", "a")
        writedlm(out0, μ)
        close(out0) 
        out1 = open(pwd()*"/betaOut$fixedRegSize", "a")
        writedlm(out1, tempBetaVec')
        close(out1)
        out2 = open(pwd()*"/varBetaOut$fixedRegSize", "a")
        writedlm(out2, varBeta')
        close(out2)
        out3 = open(pwd()*"/varEOut$fixedRegSize", "a")
        writedlm(out3, varE)
        close(out3)
        varUhat = var(X*tempBetaVec)
        out4 = open(pwd()*"/varUhatOut$fixedRegSize", "a")
        writedlm(out4, varUhat)
        close(out4)
        if onScreen==true
#            varU = var(X*tempBetaVec)
            @printf("iter %s varUhat %.2f varE %.2f\n", iter, varUhat, varE)
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end

function fileControlSt(fixedRegSize)
    for f in ["muOut$fixedRegSize" "betaOut$fixedRegSize" "varBetaOut$fixedRegSize" "varEOut$fixedRegSize" "varUhatOut$fixedRegSize"]
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function fileControl(nTraits,fixedRegSize)
    files2Remove = ["muOutMT$fixedRegSize", "varEOutMT$fixedRegSize", "covBetaOutMT$fixedRegSize", "varUOutMT$fixedRegSize"]
    for t in 1:nTraits
        push!(files2Remove,"beta"*"$t"*"Out$fixedRegSize")
        push!(files2Remove,"varBeta"*"$t"*"Out$fixedRegSize")
    end
    for f in files2Remove
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function fileControlSt2(fixedRegSize)
    for f in ["muOut$fixedRegSize" "beta1Out$fixedRegSize" "beta2Out$fixedRegSize" "beta3Out$fixedRegSize" "beta4Out$fixedRegSize" "covBetaOut$fixedRegSize" "varEOut$fixedRegSize"]
        if isfile(f)==true
            rm(f)
            println("$f removed")
        end
    end
end

function outputControl2(nRandComp,onScreen,iter,these2Keep,tempBetaMat,μ,covBeta,varE,fixedRegSize,nRegions)
    if iter in these2Keep
        out0 = open(pwd()*"/muOut$fixedRegSize", "a")
        writedlm(out0, μ)
        close(out0)
        for t in 1:nRandComp
            out1 = open(pwd()*"/beta"*"$t"*"Out$fixedRegSize", "a")
            writedlm(out1, tempBetaMat[t,:]')
            close(out1)
        end
        outCov = open(pwd()*"/covBetaOut$fixedRegSize", "a")
        printThis = [vcat(covBeta[r]...) for r in 1:nRegions]'
        writedlm(outCov, printThis)
        close(outCov)
        out3 = open(pwd()*"/varEOut$fixedRegSize", "a")
        writedlm(out3, varE)
        close(out3)  
        if onScreen==true
            coVarBeta = cov(tempBetaMat')
            corBeta   = cor(tempBetaMat') 
            println("iter $iter \n coVarBeta (Overall): $coVarBeta \n corBeta: $corBeta \n varE: $varE \n")
        elseif onScreen==false
             @printf("iter %s\n", iter)
        end
    end
end


function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
end

function sampleVarBeta(νS_β,whichLoci,df_β,regionSize)
    return((νS_β + dot(whichLoci,whichLoci))/rand(Chisq(df_β + regionSize)))
end
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end
function sampleVarE_w(νS_e,yCorVec,wVec,df_e,nRecords)
    return((νS_e + dot((yCorVec.*wVec),yCorVec))/rand(Chisq(df_e + nRecords)))
end
function sampleCovBeta(dfβ, regionSize, Vb , tempBetaMat, theseLoci)
    Sb = tempBetaMat[:,theseLoci]*tempBetaMat[:,theseLoci]'
#    println("Sb: $(Sb) \n Vb: $(Vb)")
    return rand(InverseWishart(dfβ + regionSize, Vb + Sb))
end

###Experimental iW ____ ALMOST DOUBLED THE SPEED
function sampleCovBeta_iW(dfβ, regionSize, Vb , tempBetaMat, theseLoci)
    Sb = tempBetaMat[:,theseLoci]*tempBetaMat[:,theseLoci]'
    return inv(rand(Wishart(dfβ + regionSize, convert(Array,Symmetric(inv(Vb + Sb))))))
end
