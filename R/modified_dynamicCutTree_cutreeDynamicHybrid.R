# copy of dynamicTreeCutHybrid function, adjusted to return the full results.
# https://github.com/cran/dynamicTreeCut/tree/master/R 
# cutreeHybrid --------------------------------------------------------------
cutreeHybrid_fullres <- 
function (dendro, distM, cutHeight = NULL, minClusterSize = 20, 
          deepSplit = 1, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, 
          minAbsGap = NULL, minSplitHeight = NULL, minAbsSplitHeight = NULL, 
          externalBranchSplitFnc = NULL, minExternalSplit = NULL, externalSplitOptions = list(), 
          externalSplitFncNeedsDistance = NULL, assumeSimpleExternalSpecification = TRUE, 
          pamStage = TRUE, pamRespectsDendro = TRUE, useMedoids = FALSE, 
          maxPamDist = cutHeight, respectSmallClusters = TRUE, verbose = 2, 
          indent = 0) 
{
  spaces = indentSpaces(indent)
  nMerge = length(dendro$height)
  if (nMerge < 1) 
    stop("The given dendrogram is suspicious: number of merges is zero.")
  if (is.null(distM)) 
    stop("distM must be non-NULL")
  if (is.null(dim(distM))) 
    stop("distM must be a matrix.")
  if ((dim(distM)[1] != nMerge + 1) | (dim(distM)[2] != nMerge + 
                                       1)) 
    stop("distM has incorrect dimensions.")
  if (pamRespectsDendro & !respectSmallClusters) 
    printFlush(paste("cutreeHybrid Warning: parameters pamRespectsDendro (TRUE)", 
                     "and respectSmallClusters (FALSE) imply contradictory intent.\n", 
                     "Although the code will work, please check you really", 
                     "intented these settings for the two arguments."))
  if (any(diag(distM != 0))) 
    diag(distM) = 0
  refQuantile = 0.05
  refMerge = round(nMerge * refQuantile)
  if (refMerge < 1) 
    refMerge = 1
  refHeight = dendro$height[refMerge]
  if (is.null(cutHeight)) {
    cutHeight = 0.99 * (max(dendro$height) - refHeight) + 
      refHeight
    if (verbose > 0) 
      printFlush(paste(spaces, "..cutHeight not given, setting it to", 
                       signif(cutHeight, 3), " ===>  99% of the (truncated) height range in dendro."))
  }
  else {
    if (cutHeight > max(dendro$height)) 
      cutHeight = max(dendro$height)
  }
  if (is.null(maxPamDist)) 
    maxPamDist = cutHeight
  nMergeBelowCut = sum(dendro$height <= cutHeight)
  if (nMergeBelowCut < minClusterSize) {
    if (verbose > 0) 
      printFlush(paste(spaces, "cutHeight set too low: no merges below the cut."))
    return(list(labels = rep(0, times = nMerge + 1)))
  }
  nExternalSplits = length(externalBranchSplitFnc)
  if (nExternalSplits > 0) {
    if (length(minExternalSplit) < 1) 
      stop("'minExternalBranchSplit' must be given.")
    if (assumeSimpleExternalSpecification && nExternalSplits == 
        1) {
      externalSplitOptions = list(externalSplitOptions)
    }
    externalBranchSplitFnc = lapply(externalBranchSplitFnc, 
                                    match.fun)
    for (es in 1:nExternalSplits) {
      externalSplitOptions[[es]]$tree = dendro
      if (length(externalSplitFncNeedsDistance) == 0 || 
          externalSplitFncNeedsDistance[es]) 
        externalSplitOptions[[es]]$dissimMat = distM
    }
  }
  MxBranches = nMergeBelowCut
  branch.isBasic = rep(TRUE, MxBranches)
  branch.isTopBasic = rep(TRUE, MxBranches)
  branch.failSize = rep(FALSE, MxBranches)
  branch.rootHeight = rep(NA, MxBranches)
  branch.size = rep(2, MxBranches)
  branch.nMerge = rep(1, MxBranches)
  branch.nSingletons = rep(2, MxBranches)
  branch.nBasicClusters = rep(0, MxBranches)
  branch.mergedInto = rep(0, MxBranches)
  branch.attachHeight = rep(NA, MxBranches)
  branch.singletons = vector(mode = "list", length = MxBranches)
  branch.basicClusters = vector(mode = "list", length = MxBranches)
  branch.mergingHeights = vector(mode = "list", length = MxBranches)
  branch.singletonHeights = vector(mode = "list", length = MxBranches)
  nBranches = 0
  spyIndex = NULL
  if (file.exists(".dynamicTreeCutSpyFile")) {
    spyIndex = read.table(".dynamicTreeCutSpyFile", 
                          header = FALSE)
    printFlush("Found 'spy file' with indices of objects to watch for.")
    spyIndex = as.numeric(spyIndex[, 1])
    printFlush(paste(spyIndex, collapse = ", "))
  }
  defMCS = c(0.64, 0.73, 0.82, 0.91, 0.95)
  defMG = (1 - defMCS) * 3/4
  nSplitDefaults = length(defMCS)
  if (is.logical(deepSplit)) 
    deepSplit = as.integer(deepSplit) * (nSplitDefaults - 
                                           2)
  deepSplit = deepSplit + 1
  if ((deepSplit < 1) | (deepSplit > nSplitDefaults)) 
    stop(paste("Parameter deepSplit (value", deepSplit, 
               ") out of range: allowable range is 0 through", 
               nSplitDefaults - 1))
  if (is.null(maxCoreScatter)) 
    maxCoreScatter = .interpolate(defMCS, deepSplit)
  if (is.null(minGap)) 
    minGap = .interpolate(defMG, deepSplit)
  if (is.null(maxAbsCoreScatter)) 
    maxAbsCoreScatter = refHeight + maxCoreScatter * (cutHeight - 
                                                        refHeight)
  if (is.null(minAbsGap)) 
    minAbsGap = minGap * (cutHeight - refHeight)
  if (is.null(minSplitHeight)) 
    minSplitHeight = 0
  if (is.null(minAbsSplitHeight)) 
    minAbsSplitHeight = refHeight + minSplitHeight * (cutHeight - 
                                                        refHeight)
  nPoints = nMerge + 1
  IndMergeToBranch = rep(0, times = nMerge)
  onBranch = rep(0, nPoints)
  RootBranch = 0
  if (verbose > 2) {
    printFlush(paste(spaces, "..Going through the merge tree"))
    pind = .initProgInd()
  }
  mergeDiagnostics = data.frame(smI = rep(NA, nMerge), smSize = rep(NA, nMerge), smCrSc = rep(NA, nMerge), smGap = rep(NA, nMerge), 
                                lgI = rep(NA, nMerge), lgSize = rep(NA, nMerge), lgCrSc = rep(NA, nMerge), lgGap = rep(NA, nMerge), 
                                merged = rep(NA, nMerge))
  if (nExternalSplits > 0) {
    externalMergeDiags = matrix(NA, nMerge, nExternalSplits)
    colnames(externalMergeDiags) = paste("externalBranchSplit", 
                                         1:nExternalSplits, sep = ".")
  }
  extender = rep(0, .chunkSize)
  for (merge in 1:nMerge) if (dendro$height[merge] <= cutHeight) {
    if (dendro$merge[merge, 1] < 0 & dendro$merge[merge, 
                                                  2] < 0) {
      nBranches = nBranches + 1
      branch.isBasic[nBranches] = TRUE
      branch.isTopBasic[nBranches] = TRUE
      branch.singletons[[nBranches]] = c(-dendro$merge[merge, 
                                                       ], extender)
      branch.basicClusters[[nBranches]] = extender
      branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                                 2), extender)
      branch.singletonHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                                   2), extender)
      IndMergeToBranch[merge] = nBranches
      RootBranch = nBranches
    }
    else if (sign(dendro$merge[merge, 1]) * sign(dendro$merge[merge, 
                                                              2]) < 0) {
      clust = IndMergeToBranch[max(dendro$merge[merge, 
                                                ])]
      if (clust == 0) 
        stop("Internal error: a previous merge has no associated cluster. Sorry!")
      gene = -min(dendro$merge[merge, ])
      ns = branch.nSingletons[clust] + 1
      nm = branch.nMerge[clust] + 1
      if (branch.isBasic[clust]) {
        if (ns > length(branch.singletons[[clust]])) {
          branch.singletons[[clust]] = c(branch.singletons[[clust]], 
                                         extender)
          branch.singletonHeights[[clust]] = c(branch.singletonHeights[[clust]], 
                                               extender)
        }
        branch.singletons[[clust]][ns] = gene
        branch.singletonHeights[[clust]][ns] = dendro$height[merge]
      }
      else {
        onBranch[gene] = clust
      }
      if (nm >= length(branch.mergingHeights[[clust]])) 
        branch.mergingHeights[[clust]] = c(branch.mergingHeights[[clust]], 
                                           extender)
      branch.mergingHeights[[clust]][nm] = dendro$height[merge]
      branch.size[clust] = branch.size[clust] + 1
      branch.nMerge[clust] = nm
      branch.nSingletons[clust] = ns
      IndMergeToBranch[merge] = clust
      RootBranch = clust
    }
    else {
      clusts = IndMergeToBranch[dendro$merge[merge, ]]
      sizes = branch.size[clusts]
      rnk = rank(sizes, ties.method = "first")
      small = clusts[rnk[1]]
      large = clusts[rnk[2]]
      sizes = sizes[rnk]
      branch1 = branch.singletons[[large]][1:sizes[2]]
      branch2 = branch.singletons[[small]][1:sizes[1]]
      spyMatch = FALSE
      if (!is.null(spyIndex)) {
        n1 = length(intersect(branch1, spyIndex))
        if ((n1/length(branch1) > 0.99 && n1/length(spyIndex) > 
             0.99)) {
          printFlush(paste("Found spy match for branch 1 on merge", 
                           merge))
          spyMatch = TRUE
        }
        n2 = length(intersect(branch2, spyIndex))
        if ((n2/length(branch1) > 0.99 && n2/length(spyIndex) > 
             0.99)) {
          printFlush(paste("Found spy match for branch 2 on merge", 
                           merge))
          spyMatch = TRUE
        }
      }
      if (branch.isBasic[small]) {
        coresize = .CoreSize(branch.nSingletons[small], 
                             minClusterSize)
        Core = branch.singletons[[small]][c(1:coresize)]
        SmAveDist = mean(colSums(distM[Core, Core, drop = FALSE])/(coresize - 
                                                                     1))
      }
      else {
        SmAveDist = 0
      }
      if (branch.isBasic[large]) {
        coresize = .CoreSize(branch.nSingletons[large], 
                             minClusterSize)
        Core = branch.singletons[[large]][c(1:coresize)]
        LgAveDist = mean(colSums(distM[Core, Core])/(coresize - 
                                                       1))
      }
      else {
        LgAveDist = 0
      }
      mergeDiagnostics[merge, ] = c(small, branch.size[small], 
                                    SmAveDist, dendro$height[merge] - SmAveDist, 
                                    large, branch.size[large], LgAveDist, dendro$height[merge] - 
                                      LgAveDist, NA)
      SmallerScores = c(branch.isBasic[small], branch.size[small] < 
                          minClusterSize, SmAveDist > maxAbsCoreScatter, 
                        dendro$height[merge] - SmAveDist < minAbsGap, 
                        dendro$height[merge] < minAbsSplitHeight)
      if (SmallerScores[1] * sum(SmallerScores[-1]) > 0) {
        DoMerge = TRUE
        SmallerFailSize = !(SmallerScores[3] | SmallerScores[4])
      }
      else {
        LargerScores = c(branch.isBasic[large], branch.size[large] < 
                           minClusterSize, LgAveDist > maxAbsCoreScatter, 
                         dendro$height[merge] - LgAveDist < minAbsGap, 
                         dendro$height[merge] < minAbsSplitHeight)
        if (LargerScores[1] * sum(LargerScores[-1]) > 
            0) {
          DoMerge = TRUE
          SmallerFailSize = !(LargerScores[3] | LargerScores[4])
          x = small
          small = large
          large = x
          sizes = rev(sizes)
        }
        else {
          DoMerge = FALSE
        }
      }
      if (DoMerge) {
        mergeDiagnostics$merged[merge] = 1
      }
      if (!DoMerge && (nExternalSplits > 0) && branch.isBasic[small] && 
          branch.isBasic[large]) {
        if (verbose > 4) 
          printFlush(paste0("Entering external split code on merge ", 
                            merge))
        branch1 = branch.singletons[[large]][1:sizes[2]]
        branch2 = branch.singletons[[small]][1:sizes[1]]
        if (verbose > 4 | spyMatch) 
          printFlush(paste0("  ..branch lengths: ", 
                            sizes[1], ", ", sizes[2]))
        es = 0
        while (es < nExternalSplits && !DoMerge) {
          es = es + 1
          args = externalSplitOptions[[es]]
          args = c(args, list(branch1 = branch1, branch2 = branch2))
          extSplit = do.call(externalBranchSplitFnc[[es]], 
                             args)
          if (spyMatch) 
            printFlush(" .. external criterion ", 
                       es, ": ", extSplit)
          DoMerge = extSplit < minExternalSplit[es]
          externalMergeDiags[merge, es] = extSplit
          mergeDiagnostics$merged[merge] = if (DoMerge) 
            2
          else 0
        }
      }
      if (DoMerge) {
        branch.failSize[[small]] = SmallerFailSize
        branch.mergedInto[small] = large
        branch.attachHeight[small] = dendro$height[merge]
        branch.isTopBasic[small] = FALSE
        nss = branch.nSingletons[small]
        nsl = branch.nSingletons[large]
        ns = nss + nsl
        if (branch.isBasic[large]) {
          nExt = ceiling((ns - length(branch.singletons[[large]]))/.chunkSize)
          if (nExt > 0) {
            if (verbose > 5) 
              printFlush(paste("Extending singletons for branch", 
                               large, "by", nExt, " extenders."))
            branch.singletons[[large]] = c(branch.singletons[[large]], 
                                           rep(extender, nExt))
            branch.singletonHeights[[large]] = c(branch.singletonHeights[[large]], 
                                                 rep(extender, nExt))
          }
          branch.singletons[[large]][(nsl + 1):ns] = branch.singletons[[small]][1:nss]
          branch.singletonHeights[[large]][(nsl + 1):ns] = branch.singletonHeights[[small]][1:nss]
          branch.nSingletons[large] = ns
        }
        else {
          if (!branch.isBasic[small]) 
            stop("Internal error: merging two composite clusters. Sorry!")
          onBranch[branch.singletons[[small]]] = large
        }
        nm = branch.nMerge[large] + 1
        if (nm > length(branch.mergingHeights[[large]])) 
          branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], 
                                             extender)
        branch.mergingHeights[[large]][nm] = dendro$height[merge]
        branch.nMerge[large] = nm
        branch.size[large] = branch.size[small] + branch.size[large]
        IndMergeToBranch[merge] = large
        RootBranch = large
      }
      else {
        if (branch.isBasic[large] & !branch.isBasic[small]) {
          x = large
          large = small
          small = x
          sizes = rev(sizes)
        }
        if (branch.isBasic[large] | (pamStage & pamRespectsDendro)) {
          nBranches = nBranches + 1
          branch.attachHeight[c(large, small)] = dendro$height[merge]
          branch.mergedInto[c(large, small)] = nBranches
          if (branch.isBasic[small]) {
            addBasicClusters = small
          }
          else addBasicClusters = branch.basicClusters[[small]]
          if (branch.isBasic[large]) {
            addBasicClusters = c(addBasicClusters, large)
          }
          else addBasicClusters = c(addBasicClusters, 
                                    branch.basicClusters[[large]])
          branch.isBasic[nBranches] = FALSE
          branch.isTopBasic[nBranches] = FALSE
          branch.basicClusters[[nBranches]] = addBasicClusters
          branch.mergingHeights[[nBranches]] = c(rep(dendro$height[merge], 
                                                     2), extender)
          branch.nMerge[nBranches] = 2
          branch.size[nBranches] = sum(sizes)
          branch.nBasicClusters[nBranches] = length(addBasicClusters)
          IndMergeToBranch[merge] = nBranches
          RootBranch = nBranches
        }
        else {
          addBasicClusters = if (branch.isBasic[small]) 
            small
          else branch.basicClusters[[small]]
          nbl = branch.nBasicClusters[large]
          nb = branch.nBasicClusters[large] + length(addBasicClusters)
          if (nb > length(branch.basicClusters[[large]])) {
            nExt = ceiling((nb - length(branch.basicClusters[[large]]))/.chunkSize)
            branch.basicClusters[[large]] = c(branch.basicClusters[[large]], 
                                              rep(extender, nExt))
          }
          branch.basicClusters[[large]][(nbl + 1):nb] = addBasicClusters
          branch.nBasicClusters[large] = nb
          branch.size[large] = branch.size[large] + branch.size[small]
          nm = branch.nMerge[large] + 1
          if (nm > length(branch.mergingHeights[[large]])) 
            branch.mergingHeights[[large]] = c(branch.mergingHeights[[large]], 
                                               extender)
          branch.mergingHeights[[large]][nm] = dendro$height[merge]
          branch.nMerge[large] = nm
          branch.attachHeight[small] = dendro$height[merge]
          branch.mergedInto[small] = large
          IndMergeToBranch[merge] = large
          RootBranch = large
        }
      }
    }
    if (verbose > 2) 
      pind = .updateProgInd(merge/nMerge, pind)
  }
  if (verbose > 2) {
    pind = .updateProgInd(1, pind)
    printFlush("")
  }
  if (verbose > 2) 
    printFlush(paste(spaces, "..Going through detected branches and marking clusters.."))
  isCluster = rep(FALSE, times = nBranches)
  SmallLabels = rep(0, times = nPoints)
  for (clust in 1:nBranches) {
    if (is.na(branch.attachHeight[clust])) 
      branch.attachHeight[clust] = cutHeight
    if (branch.isTopBasic[clust]) {
      coresize = .CoreSize(branch.nSingletons[clust], minClusterSize)
      Core = branch.singletons[[clust]][c(1:coresize)]
      CoreScatter = mean(colSums(distM[Core, Core])/(coresize - 
                                                       1))
      isCluster[clust] = branch.isTopBasic[clust] & (branch.size[clust] >= 
                                                       minClusterSize) & (CoreScatter < maxAbsCoreScatter) & 
        (branch.attachHeight[clust] - CoreScatter > minAbsGap)
    }
    else {
      CoreScatter = 0
    }
    if (branch.failSize[clust]) 
      SmallLabels[branch.singletons[[clust]]] = clust
  }
  if (!respectSmallClusters) 
    SmallLabels = rep(0, times = nPoints)
  if (verbose > 2) 
    printFlush(paste(spaces, "..Assigning Tree Cut stage labels.."))
  Colors = rep(0, times = nPoints)
  coreLabels = rep(0, times = nPoints)
  clusterBranches = c(1:nBranches)[isCluster]
  branchLabels = rep(0, nBranches)
  color = 0
  for (clust in clusterBranches) {
    color = color + 1
    Colors[branch.singletons[[clust]]] = color
    SmallLabels[branch.singletons[[clust]]] = 0
    coresize = .CoreSize(branch.nSingletons[clust], minClusterSize)
    Core = branch.singletons[[clust]][c(1:coresize)]
    coreLabels[Core] = color
    branchLabels[clust] = color
  }
  Labeled = c(1:nPoints)[Colors != 0]
  Unlabeled = c(1:nPoints)[Colors == 0]
  nUnlabeled = length(Unlabeled)
  UnlabeledExist = (nUnlabeled > 0)
  if (length(Labeled) > 0) {
    LabelFac = factor(Colors[Labeled])
    nProperLabels = nlevels(LabelFac)
  }
  else nProperLabels = 0
  if (pamStage & UnlabeledExist & nProperLabels > 0) {
    if (verbose > 2) 
      printFlush(paste(spaces, "..Assigning PAM stage labels.."))
    nPAMed = 0
    if (useMedoids) {
      Medoids = rep(0, times = nProperLabels)
      ClusterRadii = rep(0, times = nProperLabels)
      for (cluster in 1:nProperLabels) {
        InCluster = c(1:nPoints)[Colors == cluster]
        DistInCluster = distM[InCluster, InCluster]
        DistSums = colSums(DistInCluster)
        Medoids[cluster] = InCluster[which.min(DistSums)]
        ClusterRadii[cluster] = max(DistInCluster[, which.min(DistSums)])
      }
      if (respectSmallClusters) {
        FSmallLabels = factor(SmallLabels)
        SmallLabLevs = as.numeric(levels(FSmallLabels))
        nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1] == 
                                                    0)
        if (nSmallClusters > 0) 
          for (sclust in SmallLabLevs[SmallLabLevs != 
                                      0]) {
            InCluster = c(1:nPoints)[SmallLabels == sclust]
            if (pamRespectsDendro) {
              onBr = unique(onBranch[InCluster])
              if (length(onBr) > 1) 
                stop(paste("Internal error: objects in a small cluster are marked to belong", 
                           "\nto several large branches:", 
                           paste(onBr, collapse = ", ")))
              if (onBr > 0) {
                basicOnBranch = branch.basicClusters[[onBr]]
                labelsOnBranch = branchLabels[basicOnBranch]
              }
              else {
                labelsOnBranch = NULL
              }
            }
            else {
              labelsOnBranch = c(1:nProperLabels)
            }
            DistInCluster = distM[InCluster, InCluster, 
                                  drop = FALSE]
            if (length(labelsOnBranch) > 0) {
              if (length(InCluster) > 1) {
                DistSums = apply(DistInCluster, 2, sum)
                smed = InCluster[which.min(DistSums)]
                DistToMeds = distM[Medoids[labelsOnBranch], 
                                   smed]
                closest = which.min(DistToMeds)
                DistToClosest = DistToMeds[closest]
                closestLabel = labelsOnBranch[closest]
                if ((DistToClosest < ClusterRadii[closestLabel]) | 
                    (DistToClosest < maxPamDist)) {
                  Colors[InCluster] = closestLabel
                  nPAMed = nPAMed + length(InCluster)
                }
                else Colors[InCluster] = -1
              }
            }
            else Colors[InCluster] = -1
          }
      }
      Unlabeled = c(1:nPoints)[Colors == 0]
      if (length(Unlabeled > 0)) 
        for (obj in Unlabeled) {
          if (pamRespectsDendro) {
            onBr = onBranch[obj]
            if (onBr > 0) {
              basicOnBranch = branch.basicClusters[[onBr]]
              labelsOnBranch = branchLabels[basicOnBranch]
            }
            else {
              labelsOnBranch = NULL
            }
          }
          else {
            labelsOnBranch = c(1:nProperLabels)
          }
          if (!is.null(labelsOnBranch)) {
            UnassdToMedoidDist = distM[Medoids[labelsOnBranch], 
                                       obj]
            nearest = which.min(UnassdToMedoidDist)
            NearestCenterDist = UnassdToMedoidDist[nearest]
            nearestMed = labelsOnBranch[nearest]
            if ((NearestCenterDist < ClusterRadii[nearestMed]) | 
                (NearestCenterDist < maxPamDist)) {
              Colors[obj] = nearestMed
              nPAMed = nPAMed + 1
            }
          }
        }
      UnlabeledExist = (sum(Colors == 0) > 0)
    }
    else {
      ClusterDiam = rep(0, times = nProperLabels)
      for (cluster in 1:nProperLabels) {
        InCluster = c(1:nPoints)[Colors == cluster]
        nInCluster = length(InCluster)
        DistInCluster = distM[InCluster, InCluster]
        if (nInCluster > 1) {
          AveDistInClust = colSums(DistInCluster)/(nInCluster - 
                                                     1)
          ClusterDiam[cluster] = max(AveDistInClust)
        }
        else {
          ClusterDiam[cluster] = 0
        }
      }
      ColorsX = Colors
      if (respectSmallClusters) {
        FSmallLabels = factor(SmallLabels)
        SmallLabLevs = as.numeric(levels(FSmallLabels))
        nSmallClusters = nlevels(FSmallLabels) - (SmallLabLevs[1] == 
                                                    0)
        if (nSmallClusters > 0) {
          if (pamRespectsDendro) {
            for (sclust in SmallLabLevs[SmallLabLevs != 
                                        0]) {
              InCluster = c(1:nPoints)[SmallLabels == 
                                         sclust]
              onBr = unique(onBranch[InCluster])
              if (length(onBr) > 1) 
                stop(paste("Internal error: objects in a small cluster are marked to belong", 
                           "\nto several large branches:", 
                           paste(onBr, collapse = ", ")))
              if (onBr > 0) {
                basicOnBranch = branch.basicClusters[[onBr]]
                labelsOnBranch = branchLabels[basicOnBranch]
                useObjects = ColorsX %in% unique(labelsOnBranch)
                DistSClustClust = distM[InCluster, useObjects, 
                                        drop = FALSE]
                MeanDist = colMeans(DistSClustClust)
                useColorsFac = factor(ColorsX[useObjects])
                MeanMeanDist = tapply(MeanDist, useColorsFac, 
                                      mean)
                nearest = which.min(MeanMeanDist)
                NearestDist = MeanMeanDist[nearest]
                nearestLabel = as.numeric(levels(useColorsFac)[nearest])
                if (((NearestDist < ClusterDiam[nearestLabel]) | 
                     (NearestDist < maxPamDist))) {
                  Colors[InCluster] = nearestLabel
                  nPAMed = nPAMed + length(InCluster)
                }
                else Colors[InCluster] = -1
              }
            }
          }
          else {
            labelsOnBranch = c(1:nProperLabels)
            useObjects = c(1:nPoints)[ColorsX != 0]
            for (sclust in SmallLabLevs[SmallLabLevs != 
                                        0]) {
              InCluster = c(1:nPoints)[SmallLabels == 
                                         sclust]
              DistSClustClust = distM[InCluster, useObjects, 
                                      drop = FALSE]
              MeanDist = colMeans(DistSClustClust)
              useColorsFac = factor(ColorsX[useObjects])
              MeanMeanDist = tapply(MeanDist, useColorsFac, 
                                    mean)
              nearest = which.min(MeanMeanDist)
              NearestDist = MeanMeanDist[nearest]
              nearestLabel = as.numeric(levels(useColorsFac)[nearest])
              if (((NearestDist < ClusterDiam[nearestLabel]) | 
                   (NearestDist < maxPamDist))) {
                Colors[InCluster] = nearestLabel
                nPAMed = nPAMed + length(InCluster)
              }
              else Colors[InCluster] = -1
            }
          }
        }
      }
      Unlabeled = c(1:nPoints)[Colors == 0]
      if (length(Unlabeled) > 0) {
        if (pamRespectsDendro) {
          unlabOnBranch = Unlabeled[onBranch[Unlabeled] > 
                                      0]
          for (obj in unlabOnBranch) {
            onBr = onBranch[obj]
            basicOnBranch = branch.basicClusters[[onBr]]
            labelsOnBranch = branchLabels[basicOnBranch]
            useObjects = ColorsX %in% unique(labelsOnBranch)
            useColorsFac = factor(ColorsX[useObjects])
            UnassdToClustDist = tapply(distM[useObjects, 
                                             obj], useColorsFac, mean)
            nearest = which.min(UnassdToClustDist)
            NearestClusterDist = UnassdToClustDist[nearest]
            nearestLabel = as.numeric(levels(useColorsFac)[nearest])
            if ((NearestClusterDist < ClusterDiam[nearestLabel]) | 
                (NearestClusterDist < maxPamDist)) {
              Colors[obj] = nearestLabel
              nPAMed = nPAMed + 1
            }
          }
        }
        else {
          useObjects = c(1:nPoints)[ColorsX != 0]
          useColorsFac = factor(ColorsX[useObjects])
          nUseColors = nlevels(useColorsFac)
          UnassdToClustDist = apply(distM[useObjects, 
                                          Unlabeled, drop = FALSE], 2, tapply, useColorsFac, 
                                    mean)
          dim(UnassdToClustDist) = c(nUseColors, length(Unlabeled))
          nearest = apply(UnassdToClustDist, 2, which.min)
          nearestDist = apply(UnassdToClustDist, 2, min)
          nearestLabel = as.numeric(levels(useColorsFac)[nearest])
          assign = (nearestDist < ClusterDiam[nearestLabel]) | 
            (nearestDist < maxPamDist)
          Colors[Unlabeled[assign]] = nearestLabel[assign]
          nPAMed = nPAMed + sum(assign)
        }
      }
    }
    if (verbose > 2) 
      printFlush(paste(spaces, "....assigned", nPAMed, 
                       "objects to existing clusters."))
  }
  Colors[Colors < 0] = 0
  UnlabeledExist = (sum(Colors == 0) > 0)
  NumLabs = as.numeric(as.factor(Colors))
  Sizes = table(NumLabs)
  if (UnlabeledExist) {
    if (length(Sizes) > 1) {
      SizeRank = c(1, rank(-Sizes[2:length(Sizes)], ties.method = "first") + 
                     1)
    }
    else {
      SizeRank = 1
    }
    OrdNumLabs = SizeRank[NumLabs]
  }
  else {
    SizeRank = rank(-Sizes[1:length(Sizes)], ties.method = "first")
    OrdNumLabs = SizeRank[NumLabs]
  }
  ordCoreLabels = OrdNumLabs - UnlabeledExist
  ordCoreLabels[coreLabels == 0] = 0
  if (verbose > 0) 
    printFlush(paste(spaces, "..done."))
  
  return(list(labels = OrdNumLabs - UnlabeledExist, 
              cores = ordCoreLabels, 
              smallLabels = SmallLabels, 
              onBranch = onBranch, 
              mergeDiagnostics = if (nExternalSplits ==  0) mergeDiagnostics else cbind(mergeDiagnostics,  externalMergeDiags), 
              mergeCriteria = list(maxCoreScatter = maxCoreScatter, 
              minGap = minGap, maxAbsCoreScatter = maxAbsCoreScatter, 
              minAbsGap = minAbsGap, minExternalSplit = minExternalSplit), 
              branches = list(nBranches = nBranches, IndMergeToBranch = IndMergeToBranch, 
              RootBranch = RootBranch, isCluster = isCluster, nPoints = nMerge + 1)))
}

# other internal functions https://github.com/cran/dynamicTreeCut/blob/master/R/treeCut.R

# Progress indicator...

.initProgInd = function( leadStr = "..", trailStr = "", quiet = !interactive())
{
  oldStr = " ";
  cat(oldStr);
  progInd = list(oldStr = oldStr, leadStr = leadStr, trailStr = trailStr);
  class(progInd) = "progressIndicator";
  .updateProgInd(0, progInd, quiet);
}

.updateProgInd = function(newFrac, progInd, quiet = !interactive())
{
  if (class(progInd)!="progressIndicator")
    stop("Parameter progInd is not of class 'progressIndicator'. Use initProgInd() to initialize",
         "it prior to use.");
  
  newStr = paste(progInd$leadStr, as.integer(newFrac*100), "% ", progInd$trailStr, sep = "");
  if (newStr!=progInd$oldStr)
  {
    if (quiet)
    {
      progInd$oldStr = newStr;
    } else {
      cat(paste(rep("\b", nchar(progInd$oldStr)), collapse=""));
      cat(newStr);
      if (exists("flush.console")) flush.console();
      progInd$oldStr = newStr;
    }
  }
  progInd;
}


# The following are supporting function for GetClusters. 

.CoreSize = function(BranchSize, minClusterSize)
{
  BaseCoreSize = minClusterSize/2 + 1;
  if (BaseCoreSize < BranchSize)
  {
    CoreSize = as.integer(BaseCoreSize + sqrt(BranchSize - BaseCoreSize));
  } else CoreSize = BranchSize;
  CoreSize;
}

# This assumes the diagonal of the distance matrix
# is zero, BranchDist is a square matrix whose dimension is at least 2.

.CoreScatter = function(BranchDist, minClusterSize)
{
  nPoints = dim(BranchDist)[1];
  PointAverageDistances = colSums(BranchDist) / (nPoints-1);
  CoreSize = minClusterSize/2 + 1;
  if (CoreSize < nPoints)
  {
    EffCoreSize = as.integer(CoreSize + sqrt(nPoints - CoreSize));
    ord = order(PointAverageDistances);
    Core = ord[c(1:EffCoreSize)];
  } else {
    Core = c(1:nPoints);
    EffCoreSize = nPoints;
  }
  CoreAverageDistances = colSums(BranchDist[Core, Core]) / (EffCoreSize-1);
  mean(CoreAverageDistances);
}

.interpolate = function(data, index)
{
  i = round(index);
  n = length(data);
  if (i<1) return(data[1]);
  if (i>=n) return(data[n]);
  
  r = index-i;
  data[i] * (1-r) + data[i+1] * r;
}

.chunkSize = 100