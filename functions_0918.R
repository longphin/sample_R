suppressPackageStartupMessages(library(glmnet));
suppressPackageStartupMessages(library(elasticnet));
suppressPackageStartupMessages(library(parser));
Sys.setlocale('LC_ALL','C');

# ===============================================
# Getting results
# ===============================================


# input: give this function a vector of single words and the size for the phrases
# output: a vector of phrases, each phrase made of #size words
getPhrases<-function(text, size){
  # size = the number of words making up the phrases
  
  # error checking: text does not contain enough words to make a phrase of that size
  if(length(text)<size) return(NULL);
  
  invisible(sapply(1:(length(text)-size+1), function(n)	paste(text[n:(n+size-1)], collapse=' ')))
}


# input: give this function a single string
# output: a vector of all words separated, without punctuations and in lower case
separateToWords<-function(string, minWordSize=1, stopwords='', phraseSize=1){
  # minWordSize = only keep words that have more than that number of characters
  # phraseSize = the number of words that make up a phrase. Can be a sequence like 1:3 obtains phrases made of 1 words, 2 words, and 3 words
  
  # break up string into words
  onlyWords=unlist(strsplit(gsub("\n", " ", string), "[^[:alpha:]]"));
  onlyWords=tolower(onlyWords[onlyWords!='']);
  
  # remove words that are too small?
  if(minWordSize>1)
    onlyWords=onlyWords[which(nchar(onlyWords)>minWordSize)];
  
  # remove stopwords, if used
  if((length(stopwords)==1 && stopwords!='') || (length(stopwords)>1 && !is.na(stopwords[1])))
    onlyWords=onlyWords[onlyWords %in% tolower(stopwords) == FALSE];
  
  # error checking: if string only contained unwanted text, it should become empty by this step.
  if(length(onlyWords)==0) return(NULL);
  # error checking: if user only wanted single words (phraseSize=1), then return it now
  if(length(phraseSize)==1 && phraseSize<=1) return(onlyWords);
  
  # obtain phrases made of 2 or more words, if chosen
  onlyPhrases=unlist(lapply(phraseSize, function(size) getPhrases(onlyWords, size)));
  # remove stopwords from phrases
  onlyPhrases=onlyPhrases[onlyPhrases %in% tolower(stopwords) == FALSE];
  
  # if user wanted single words to be included in the output, append it to onlyPhrases. This occurs when phraseSize contains the number 1
  if(1 %in% phraseSize) return(c(onlyWords, onlyPhrases));
  
  # otherwise, user only wanted phrases
  invisible(onlyPhrases);
}


# input: give this function a list of vectors. The vectors contains single words.
# output: a matrix with words in the columns, rows are the elements of the list
createFeaturesMatrix<-function(List){
  # not used; (this would have removed words only used once)
  #  unimportantwords=names(which(table(unlist(List))==1))
  #  List=lapply(List, function(i) setdiff(i, unimportantwords))
  
  # find unique words/phrases
  words=unique(unlist(List));
  
  # preallocate a matrix. The frequencies will be added to it.
  mat=matrix(0, nrow=length(List), ncol=length(words));
  colnames(mat)=words;
  
  # gets frequency of each word per article
  li=lapply(List, function(i) table(i));
  
  # add the frequencies to the matrix
  for(i in 1:length(li))
  {
    for(word in names(li[[i]]))
    {
      mat[i, word]=li[[i]][word]
    }
  }
  
  invisible(mat)
}


# input: a matrix and the reweighting method to be used
#        method=="tfidf"||"l2"||"<anything else>"
#        stopwords is a vector of words to remove
# output: a reweighted matrix
reweightMatrix<-function(mat, method='l2', stopwords=''){
  # remove stopwords
  if((length(stopwords)==1 && stopwords!='') || (length(stopwords)>1 && !is.na(stopwords[1])))
  {
    mat=mat[, which(!colnames(mat)%in%stopwords)]
  }
  
  # l2 reweighting
  if(method=="l2")
    return(sweep(mat, 2, sqrt(colSums(mat^2)), "/"));
  
  # tfidf reweighting
  if(method=="tfidf")
    return(mat/rowSums(mat)*log(nrow(mat)/colSums(mat>0)));
  
  # only stopwords removal used
  invisible(mat)
}


# input: a vector of words and desired length
# output: a vector of words without subphrases and of the desired length
removeSubphrases<-function(txt, desiredlength=10)
{
  start=1;
  end=desiredlength;
  ans=head(txt, desiredlength);
  
  # remove subphrases from the terms chosen. Repeat until no subphrases are in the chosen terms.
  repeat
  {
    # remove subphrases
    for(i in 1:length(ans))
    {
      # if ans[i] is subphrase of any other words, remove it.
      if(length(grep(paste("(.* )?", ans[i], "( .*)?", sep=''), ans))>1)
      {
        ans[i]=""
      }
    }
    
    # remove empty strings
    ans=ans[ans!=""];
    
    # no subphrases are in the terms
    if(length(ans)==desiredlength) return(ans);
    
    # there were still are subphrases in the terms. Pick more terms to get to the desired length.
    start=end+1;
    end=end+desiredlength-length(ans);
    
    # no more subphrases
    if(end<start || end>length(txt)) return(ans);
    
    # subphrases remaining, continue to loop
    ans=c(ans, txt[start:end]);
  }
  
  invisible(ans)
}


# input: character vector. If an element is "file:<filename>", then it will
#        read the contents of <filename> to obtain words
# output: a vector of stopwords
getStopwords<-function(vec){
  # finding if any files are to be read
  locateFiles=grep("^file:", vec);
  filesToRead = gsub('^file:(./)?(.*)$', "\\2", vec[locateFiles]);
  
  # read any files to create the list of words
  stopwords_fromfile=unlist(lapply(filesToRead, function(file) readLines(file)));
  # combine the words from files and words no from files
  invisible(c(stopwords_fromfile, vec[setdiff(1:length(vec), locateFiles)]))
}


# input: features matrix, y vector, and number of terms desired
# output: a vector with the results from cooccurence
getCooc<-function(X, y, terms=10){
  # criterion = mean(columns given positive article)
  criterion=sort(colMeans(as.matrix(X[which(y==1),])), decreasing=TRUE);
  invisible(criterion[removeSubphrases(names(criterion), terms)])
}


# input: features matrix, y vector, and number of terms desired
# output: a vector with the results from cooccurence
getCorr<-function(X, y, terms=10){
  criterion=as.vector(abs(cor(X, y)));
  names(criterion)=colnames(X);
  criterion=sort(criterion, decreasing=TRUE);
  invisible(criterion[removeSubphrases(names(criterion), terms)])
}


# input: features matrix and y vector
# output: L1LR result
getGlmnet<-function(X, y, alpha=1, lambda=0.01, model="binomial"){ # model=binomial for l1lr, gaussian for lasso
  mo1=tryCatch({invisible(glmnet(X, y, alpha=1, pmax=15, family=model))},error=function(e) NULL);
  # error checking: empty model
  if(is.null(mo1)) return(c(' '=' '))
  
  #	lam=mo1$lambda[length(mo1$lambda)]
  #	mo2=tryCatch({invisible(glmnet(X, y, alpha=1, lambda=lam, family=model))}, error=function(e) NULL);#, maxit=2147483647)
  # error checking: empty model
  #	if(is.null(mo2)) return(c(' '=' '))
  
  #  limited=which(mo2$beta[,1]>(max(mo2$beta)*.1))
  #	re1 = sort(mo2$beta[names(limited),], decreasing=TRUE);
  limited=which(mo1$beta[,ncol(mo1$beta)]>(max(mo1$beta[,ncol(mo1$beta)])*.1));
  
  # get the resulting terms in order of significance
  re1 = sort(mo1$beta[names(limited),ncol(mo1$beta)], decreasing=TRUE);
  # fixing: if only 1 term is found, the name is not properly given
  if(length(re1)==1)
    names(re1)=names(limited);
  
  invisible(re1)
}


# input: a scaled features matrix (columns with mean 0, variance of 1), y vector
# output: result from SPCA
getSPCA<-function(X){
  SPCA=arrayspc(X, K=1, para=rep(0.1, 1));
  
  matchings=which(SPCA$loadings[, 1]!=0);
  result=SPCA$loadings[matchings];
  names(result)=colnames(X)[matchings];
  
  invisible(result)
}


# input: filename
# output: list<files> which contains list<articles>. list<articles> contains
#         vectors with names "title","author","review"
getDataFromFile<-function(filenames="masterfile")
{
  ans=lapply(filenames, function(filename)
  {
    # number of lines to read
    len=nlines(filename)
    
    con<-file(filename, 'r')
    articles=lapply(1:(len/4), function(i)
    {
      input<-readLines(con, n=4)
      
      input[1]=gsub(';', '\n', input[1]) # author row
      input=gsub('[[:punct:]]', ' ', input) # remove punctuations from author
      # create vector of author, title, and review
      tmp=c(author=input[1], title=input[2], review=input[3])
      return(tmp)
    })
    close(con)
    return(articles)
  })
  
  invisible(ans)
}


# input: data matrix (not features matrix)
# output: vector where the titles and reviews are concatenated
getStrings<-function(mat){
  invisible(paste(mat['title', ], mat['review', ])) # puts the title and review of each article into single strings
}


# input: filename containing all articles
# output: a data matrix where the rows contain the titles, authors, and reviews. Columns are article numbers
createDataMatrix<-function(filename){
  # get the text in lists
  tmp=getDataFromFile(filename);
  
  # get the authors, titles, and reviews to be put into a data matrix
  ss=tolower(unlist(tmp));
  a=ss[which(names(ss)=='author')]; # get all of the authors;
  b=ss[which(names(ss)=='title')]; # get all of the titles;
  c=ss[which(names(ss)=='review')]; # get all of the reviews;
  
  # put authors, titles, reviews into data matrix
  mat=rbind(author=a,title=b,review=c);
  colnames(mat)=1:ncol(mat); # colnames will be article number as obtained
  
  # remove duplicates articles. Duplication is based on if reviews are equal
  t1=gsub(" {2,}", " ", gsub(" *$", "", gsub("^ *", "", mat['review',])))
  isdup=duplicated(t1)
  keepuniques=which(isdup==FALSE)
  datamat=mat[,keepuniques]
  
  invisible(datamat)
}


# input: data matrix (NOT feature matrix), where the rows contain the titles,
#        authors, and reviews. Also, input queries for the title, author, review.
# output: y vector of 1's and -1's
createYVector<-function(mat, title='', author='', review=''){
  # preallocate vector
  y=vector(mode="logical", length=ncol(mat));
  
  # if title contains the query_title, indicate as TRUE
  if(title!='')
    y=y|grepl(tolower(title), mat['title',]);
  
  # if author contains query_author, indicate as TRUE
  if(author!='')
    y=y|grepl(tolower(author), mat['author',]);
  
  # if review contains query_review, indicate as TRUE
  if(review!='')
    y=y|grepl(tolower(review), mat['review',]);
  
  invisible(as.numeric(y)*2-1)
}


# input: a string and queries
# output: string with the query terms removed
removeQuery<-function(string, query_title, query_author, query_review){
  # remove query_title from titles
  if(query_title!='')
    string=gsub(query_title, '', string);
  # remove query_author from authors
  if(query_author!='')
    string=gsub(query_author, '', string);
  # remove query_review from reviews
  if(query_review!='')
    string=gsub(query_review, '', string);
  
  invisible(string)
}


# ===============================================
# functions that can be useful for demonstrations
# ===============================================
filename='alldata0909clean'
query_title=''
query_author='Loera'
query_review=''
weighting='l2'
stopwords='file:stopwords2.txt'
minchar=3
phraseSize=c(1,2)
positivesProportion=1
negativesMultiple=10
savefile=inputfile='alldata_0909.Rda'

# (NO LONGER USED)
# input: vector of characters
# output: single character in a format that checks if any of the words match
# example: vec=c("how", "you", "doing") -> "(how)|(you)|(doing)"
# convertToRegularExpression<-function(vec){
#invisible(paste('(', paste(tolower(vec), collapse=')|('), ')', sep=''))
#}


# =============
# demonstration
# =============

# input: a list obtained from demonstration()
# output: nothing, only prints the list in a nice format
printResults<-function(li){
  # making things a little easier to type
  a=li$cooc;
  a_names=names(a);
  b=li$corr;
  b_names=names(b);
  c=li$l1lr;
  c_names=names(c);
  d=li$lasso;
  d_names=names(d);
  e=li$SPCA;
  e_names=names(e);
  
  a_l=length(a);
  b_l=length(b);
  c_l=length(c);
  d_l=length(d);
  e_l=length(e);
  m=max(a_l,b_l,c_l,d_l);
  
  # make a matrix to be printed, containing the values of the criterion used
  mat=matrix('', ncol=5, nrow=m);
  # make a matrix to be printed, containing the words for each result
  mat_words=matrix('', ncol=5, nrow=m);
  for(i in 1:m){
    if(i<=a_l){
      mat[i, 1]=a[i]
      mat_words[i, 1]=a_names[i]
    }
    if(i<=b_l){
      mat[i, 2]=b[i]
      mat_words[i, 2]=b_names[i]
    }
    if(i<=c_l){
      mat[i, 3]=c[i]
      mat_words[i, 3]=c_names[i]
    }
    if(i<=d_l){
      mat[i, 4]=d[i]
      mat_words[i, 4]=d_names[i]
    }
    if(i<=e_l){
      mat[i, 5]=e[i]
      mat_words[i, 5]=e_names[i]
    }
  }
  
  # print the matrices neatly
  colnames(mat)=c("cooc", "corr", "l1lr", "lasso", "spca");
  colnames(mat_words)=c("cooc", "corr", "l1lr", "lasso", "spca");
  rownames(mat)=rownames(mat_words)=1:m;
  
  cat(paste(rep("=",46), collapse=''), '\n');
  cat(paste(rep(" ", 19), collapse=''), "Results\n");
  cat(paste(rep("=",46), collapse=''), '\n');
  cat("\n");
  
  #  xtab<-xtable(mat_words) # helps easily convert results to latex
  #  print(xtab)
  
  # print numbers
  print(strtrim(mat,10), quote=FALSE);
  
  cat("\n");
  
  # print words
  print(mat_words, quote=FALSE);
  cat("\n");
  invisible(NULL)
}
printResultsNoSPCA<-function(li, num){
  # making things a little easier to type
  a=li$cooc;
  a_names=names(a);
  b=li$corr;
  b_names=names(b);
  c=li$l1lr;
  c_names=names(c);
  d=li$lasso;
  d_names=names(d);
  
  a_l=length(a);
  b_l=length(b);
  c_l=length(c);
  d_l=length(d);
  m=max(a_l,b_l,c_l,d_l);
  
  # make a matrix to be printed, containing the words for each result
  mat_words=matrix('', ncol=4, nrow=m);
  for(i in 1:m){
    if(i<=a_l){
      mat_words[i, 1]=a_names[i]
    }
    if(i<=b_l){
      mat_words[i, 2]=b_names[i]
    }
    if(i<=c_l){
      mat_words[i, 3]=c_names[i]
    }
    if(i<=d_l){
      mat_words[i, 4]=d_names[i]
    }
  }
  
  # print the matrices neatly
  colnames(mat_words)=c("cooc", "corr", "l1lr", "lasso");
  rownames(mat_words)=1:m;
  
  cat(paste(rep("=",46), collapse=''), '\n');
  cat(paste(rep(" ", 19), collapse=''), "Results", " ", num, "\n");
  cat(paste(rep("=",46), collapse=''), '\n');
  cat("\n");
  
  #  xtab<-xtable(mat_words) # helps easily convert results to latex
  #  print(xtab)
  
  print(mat_words, quote=FALSE);
  cat("\n");
  invisible(NULL)
}



# input: stuff
# output: nothing, saves features and data matrix
createDataMatrices<-function(filename, phraseSize, minchar,
                             stopwords='', savefile='datamatrices.Rda'){
  stopwords=getStopwords(stopwords);
  
  # matrix. rows are title, author, review. cols are 1,.., # articles
  datamat=createDataMatrix(filename);
  
  # concatenate titles and reviews together
  string=paste(datamat['title',],datamat['review',]);
  
  # break up each article into their words/phrases
  words_list=lapply(string, function(article) separateToWords(article, phraseSize=phraseSize, minWordSize=minchar, stopwords=stopwords));
  
  # create features matrix
  X=createFeaturesMatrix(words_list);
  
  # put some useful matrices to be saved, so that this function never has to be
  # used ever again, which is probably for the better.
  if(savefile!=''){
    listofdatas=list(X=X, datamat=datamat)
    save(listofdatas, file=savefile);
  }
  
  invisible(NULL)
}



# use this to reduce the size of the features matrix. Requires a query
# default takes 100% of positives, and 10x that as negatives (randomly)
# output: the features matrix (sample of entire data)
createDataMatricesRandom<-function(filename, query_title='', query_author='',
                                   query_review='', phraseSize, minchar,
                                   stopwords='', positivesProportion=1, negativesMultiple=10,
                                   savefile=''){
  stopwords=getStopwords(stopwords);
  
  # data matrix
  datamat=createDataMatrix(filename);
  # y vector
  y=createYVector(datamat, query_title, query_author, query_review);
  
  positives=which(y==1)
  negatives=which(y==-1)
  
  # sample from the positives and negatives in order to work with a smaller matrix
  randomlypickpos=sample(positives, ceiling(length(positives)*positivesProportion))
  randomlypickneg=sample(negatives,
                         size=min(max(length(randomlypickpos)*negativesMultiple, 500), length(negatives)))
  
  # concatenate titles and reviews
  string=paste(datamat['title', c(randomlypickpos,randomlypickneg)]
               ,datamat['review', c(randomlypickpos,randomlypickneg)]);
  
  # break up titles+reviews into words/phrases
  words_list=lapply(string, function(article) separateToWords(article, phraseSize=phraseSize, minWordSize=minchar, stopwords=stopwords));
  
  # create feature matrix using the smaller sample
  X=createFeaturesMatrix(words_list);
  
  # save the matrices if needed
  if(savefile!=''){
    listofdatas=list(X=X, datamat=datamat[,c(randomlypickpos,randomlypickneg)])
    save(listofdatas=listofdatas, file=savefile)
  }
  
  invisible(X)
}


# input: to make demonstrations faster, the features matrix is saved in a .Rda file
#        By loading that file, we can access the already-created features matrix
# output: results
demonstration<-function(inputfile='matricesFromData.Rda', query_title='', query_author='', query_review='', weighting='l2', lambda=0.01, stopwords=''){
  load(file=inputfile)
  
  # y vector, features matrix (loaded from Rda file), and reweight the features matrix
  y=createYVector(listofdatas$datamat, query_title, query_author, query_review);
  X=listofdatas$X;
  X_rw=reweightMatrix(X, method=weighting, stopwords=stopwords);
  
  # l1lr and lasso
  l1lr = getGlmnet(X, y,lambda=lambda, model="binomial");
  lasso = getGlmnet(X, y, lambda=lambda, model="gaussian");
  # sort l1lr and lasso betas in order
  l1lr_result = sort(l1lr, decreasing=TRUE);
  lasso_result = sort(lasso, decreasing=TRUE);
  
  # spca using only positive articles
  onlyPositives=X_rw[y==1, ];
  # remove columns of all 0's, then scale the matrix
  onlyPositives=scale(onlyPositives[, which(colSums(as.matrix(onlyPositives))>0)]);
  spca_result = sort(getSPCA(onlyPositives), decreasing=TRUE);
  
  # determine how long results for coocurrence and correlation screening should
  # give. This will be the max length of l1lr/lasso
  tt=max(length(l1lr_result), length(lasso_result));
  cooc = getCooc(X_rw, y, terms=tt);
  corr = getCorr(X_rw, y, terms=tt);
  
  # make a list of results
  li=list(cooc=cooc, corr=corr, l1lr=l1lr_result, lasso=lasso_result, SPCA=spca_result, y=y);
  # print results in a nice format
  printResults(li);
  
  options(warn=0);
  invisible(li)
}


doNotUseThis<-function(){
  load('alldata_0914.Rda')
  Y=createYVector(listofdatas$datamat, query_title, query_author, query_review);
  maxlen=ncol(listofdatas$datamat)
  
  positives=which(Y==1)
  negatives=which(Y==-1)
  set.seed(12) # change this (this indicates the _r13 or _r12 at the end of result files
  randomlypickneg=sample(negatives, size=length(negatives))
  
  set.seed(13)
  randomlypickpos=sample(positives, size=floor(length(positives)/1))
  len=length(randomlypickpos)
  
  for(num in 1:1600){ # 1600 is an arbitrary number
    
    if(len*num>maxlen) break;
    seq=c(randomlypickpos,randomlypickneg[1:(len*num)])
    
    y=createYVector(listofdatas$datamat[, seq], query_title, query_author, query_review);
    X=listofdatas$X[seq, ];
    X_rw=reweightMatrix(X, method=weighting, stopwords=stopwords);
    
    l1lr = getGlmnet(X, y,lambda=lambda, model="binomial");
    lasso = getGlmnet(X, y, lambda=lambda, model="gaussian");
    
    # sort l1lr and lasso betas in order
    l1lr_result = sort(l1lr, decreasing=TRUE);
    lasso_result = sort(lasso, decreasing=TRUE);
    
    tt=max(length(l1lr_result), length(lasso_result));
    if(tt==0){ sprintf("No Result: %d",num); next;}
    cooc = getCooc(X_rw, y, terms=tt);
    corr = getCorr(X_rw, y, terms=tt);
    
    # make a list of results
    li=list(cooc=cooc, corr=corr, l1lr=l1lr_result, lasso=lasso_result);
    # print results in a nice format
    printResultsNoSPCA(li, length(seq));
  }
}
