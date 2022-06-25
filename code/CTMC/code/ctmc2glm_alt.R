ctmc2glm_alt <-
  function(ctmc, raster.list, 
           normalize.gradients=FALSE, grad.point.decreasing=FALSE,
           directions=4, zero.idx=integer()){
    
    ## Function to take a CTMC path and covariate rasters and return data
    ##  that can be fit as a Poisson GLM
    ##
    ## Inputs:
    ##
    ##  ctmc - a "ctmc" object (output from "path2ctmc.r")
    ##  stack.static - a raster stack or raster layer of "location-based" covariates
    ##  stack.grad - a raster stack or raster layer of "gradient based covariates
    ##  grad.list - A list of raster stacks, where each stack is the vector gradiant
    ##        of a covariate to be included, e.g., list(var1=var1_stack) and 
    ##        names(var1_stack) = c("grad.x","grad.y")
    ##  crw - logical.  If TRUE, then an autocovariate is computed defined as a
    ##        directional covariate pointed in the direction of the most recent
    ##        transition in the CTMC
    ##  normalize.gradients - logical.  If TRUE, then normalize all gradient
    ##        covariates by dividing by the length of the gradient vector at each point
    ##  grad.point.decreasing - logical.  If TRUE, then the gradient covariates are positive
    ##        in the direction of decreasing values of the covariate.
    ##
    examplerast = raster.list[[1]][[1]]
    
    ### Covert static rasters to gradients
    if(!is.null(raster.list$stack.grad)){
      if(!compareRaster(examplerast, raster.list$stack.grad, stopiffalse=F)) stop("All raster objects need to have the same topology!")
      stack.gradient=rast.grad(raster.list$stack.grad)
      if(normalize.gradients){
        lengths=sqrt(stack.gradient$grad.x^2+stack.gradient$grad.y^2)
        stack.gradient$grad.x <- stack.gradient$grad.x/lengths
        stack.gradient$grad.y <- stack.gradient$grad.y/lengths
      }
    } else{
      stack.gradient=list(grad.x=NULL, grad.y=NULL)
    }
    
    ### Obtain values for gradient rasters
    if(!is.null(raster.list$grad.list)){
      X.grad.x=NULL
      X.grad.y=NULL
      for(k in 1:length(raster.list$grad.list)){
        if(!compareRaster(examplerast, raster.list$grad.list[[k]], stopiffalse=F)) stop("All raster objects need to have the same topology!")
        grad.x=raster::values(raster.list$grad.list[[k]]$grad.x)
        grad.x[is.na(grad.x)] <- 0
        grad.y=raster::values(raster.list$grad.list[[k]]$grad.y)
        grad.y[is.na(grad.y)] <- 0
        X.grad.x=cbind(X.grad.x,grad.x)
        X.grad.y=cbind(X.grad.y,grad.y)  
      }
      colnames(X.grad.x) = names(raster.list$grad.list)
      colnames(X.grad.y) = names(raster.list$grad.list)
      stack.gradient2 = list(grad.x=X.grad.x, grad.y=X.grad.y)
      if(normalize.gradients){
        lengths=sqrt(stack.gradient2$grad.x^2+stack.gradient2$grad.y^2)
        stack.gradient2$grad.x <- stack.gradient2$grad.x/lengths
        stack.gradient2$grad.y <- stack.gradient2$grad.y/lengths
      }
      
    } else{
      stack.gradient2=list(grad.x=NULL, grad.y=NULL)
    }
    
    stack.gradient = list(
      grad.x = cbind(stack.gradient$grad.x, stack.gradient2$grad.x),
      grad.y = cbind(stack.gradient$grad.y, stack.gradient2$grad.y)
    )
    
    
    locs=ctmc$ec
    wait.times=ctmc$rt
    
    ##
    ## Make X matrix 
    ##
    
    ## raster cells that are NOT in "zero.idx"
    notzero.idx=1:ncell(examplerast)
    if(length(zero.idx)>0){
      notzero.idx=notzero.idx[-zero.idx]
    }
    
    ## sort.idx=sort(locs,index.return=TRUE)$ix
    ## This is for a rook's neighborhood
    ##  n.nbrs=4
    ## sort.idx=rep(sort.idx,each=n.nbrs)
    adj=adjacent(examplerast,locs,pairs=TRUE,sorted=TRUE,id=TRUE,directions=directions,target=notzero.idx)
    adj.cells=adj[,3]
    rr=rle(adj[,1])
    time.idx=rep(rr$values,times=rr$lengths)
    start.cells=adj[,2]
    z=rep(0,length(start.cells))
    idx.move=rep(0,length(z))
    diag.move=rep(0,length(locs))
    for(i in 1:(length(locs))){
      idx.t=which(time.idx==i)
      idx.m=which(adj.cells[idx.t]==locs[i+1])
      z[idx.t[idx.m]] <- 1
      if(length(idx.m)==0){
        diag.move[i]=1
      }
    }
    
    #browser()
    
    
    ## Tau
    tau=rep(wait.times,times=rr$lengths)
    ##
    t=rep(ctmc$trans.times,times=rr$lengths)
    
    ##
    ## Get x values for static covariates
    ##
    
    if(!is.null(raster.list$stack.static)){
      X.static=matrix(raster::values(raster.list$stack.static),ncol=nlayers(raster.list$stack.static))[start.cells,,drop=F]
      colnames(X.static) <- names(raster.list$stack.static)
    } else {X.static=NULL}
    
    ##
    ## Get x values for gradiant covariates
    ##
    
    xy.cell=xyFromCell(examplerast,start.cells)
    xy.adj=xyFromCell(examplerast,adj.cells)
    ## Find normalized vectors to adjacent cells
    v.adj=(xy.adj-xy.cell)/sqrt(apply((xy.cell-xy.adj)^2,1,sum))
    ## dot product for gradient covariates
    if(!(is.null(raster.list$stack.grad) & is.null(raster.list$grad.list))){
      X.grad=v.adj[,1]*stack.gradient$grad.x[start.cells,,drop=F]+v.adj[,2]*stack.gradient$grad.y[start.cells,,drop=F]
      ## Make gradient vectors point toward LOWER values (if specified)
      if(grad.point.decreasing==TRUE){
        X.grad=-X.grad
      }
      colnames(X.grad) <- paste0(colnames(stack.gradient$grad.x), "_grad")
    } else {X.grad=NULL}
    
    
    ##
    ## Get crw covariate
    ##
    
    #browser()
    
    ##
    ## should probably figure a better way to handle diagonals.
    ##
    
    ## update 20160519
    ##     idx.move=rep(,length(locs))
    ## idx.move[-which(diag.move==1)] <- which(z==1)
    idx.move=which(z==1)
    ## last point
    idx.move=c(idx.move,length(z))
    
    
    v.moves=v.adj[rep(idx.move[1:(length(rr$lengths)-1)],times=rr$lengths[-1]),]
    ## shift v.moves to be a lag 1 direction
    v.moves=rbind(matrix(0,ncol=2,nrow=rr$lengths[1]),v.moves)
    ## dot product with vectors to adjacent cells
    #browser()
    X = cbind(crw=apply(v.moves*v.adj,1,sum), north=v.adj[,2], east=v.adj[,1])
    
    #browser()
    
    X = cbind(X, X.static, X.grad)
    
    xys=cbind(xy.cell,xy.adj)
    colnames(xys)=c("x.current","y.current","x.adj","y.adj")
    X=cbind(X,xys)
    
    
    T=length(wait.times)
    p=ncol(X)
    
    #browser()
    
    
    
    out=data.frame(z=z,X,tau=tau,t=t)
    ## remove last time step
    T=nrow(out)
    out=out[-((T-3):T),]
    out
  }