pkgs <- c("quanteda","stringi","stringr","dplyr","tidyr","purrr","tibble","data.table","readr")
to_install <- pkgs[! pkgs %in% installed.packages()[,"Package"]]
if(length(to_install)>0) install.packages(to_install, dependencies=TRUE)

suppressPackageStartupMessages({
  library(quanteda)
  library(stringi)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(data.table)
  library(readr)
})

anatomy_terms <- c(
  "lung","lungs","upper lobe","middle lobe","lower lobe","right lung","left lung",
  "apex","base","hilum","perihilar","peripheral",
  "pleura","pleural space","pleural cavity","costophrenic angle","hemithorax",
  "heart","cardiac silhouette","mediastinum","aorta","trachea","carina",
  "diaphragm","subdiaphragmatic","gastric bubble",
  "rib","ribs","clavicle","scapula","sternum","spine","vertebra",
  "endotracheal tube","et tube","nasogastric tube","ng tube","central line","cvc",
  "picc line","picc","pacemaker","chest drain","catheter","electrode lead","lead"
)

observation_terms <- c(
  "opacity","opacities","consolidation","atelectasis","atelectatic change",
  "interstitial markings","reticular pattern","nodular pattern","ground glass opacity",
  "ground glass opacities","airspace disease","air space disease",
  "pleural effusion","effusion","edema","oedema","pneumothorax","pneumomediastinum",
  "subcutaneous emphysema","pneumonia","aspiration","infiltrate",
  "cardiomegaly","pulmonary congestion","vascular redistribution",
  "fracture","lytic lesion","sclerotic lesion","mass","nodule","calcification",
  "malpositioned","kinked","migrated","overlies","high position","low position",
  "hernia","scoliosis","dextrocardia"
)


context_triggers <- list(
  pre_neg = c("no","without","free of","absent","negative for","not seen","lack of"),
  post_neg = c("ruled out","excluded"),
  pre_uncert = c("possible","possibly","likely","cannot exclude","suspicious for",
                 "consistent with","in keeping with","questionable","equivocal"),
  post_uncert = c("may represent","may reflect","could represent","could reflect",
                  "suggesting","suggests","suggestive of","compatible with"),
  terminate = c("but","however","yet","although","except","though"),
  pseudo_neg = c("no change","no interval change","no significant change",
                 "no appreciable change","no substantial change")
)


rx_from_phrases <- function(phrases){
  esc <- str_replace_all(phrases,"([\\^$.|?*+()\\[\\]{}\\\\])","\\\\\\1")
  esc <- str_replace_all(esc,"\\s+","\\\\s+")
  paste0("(?i)\\b(",paste(esc,collapse="|"),")\\b")
}

rx_ctx <- list(
  pre_neg=rx_from_phrases(context_triggers$pre_neg),
  post_neg=rx_from_phrases(context_triggers$post_neg),
  pre_unc=rx_from_phrases(context_triggers$pre_uncert),
  post_unc=rx_from_phrases(context_triggers$post_uncert),
  term=rx_from_phrases(context_triggers$terminate),
  pseudo=rx_from_phrases(context_triggers$pseudo_neg)
)

locate_all <- function(text, rx){
  m <- stringi::stri_locate_all_regex(text, rx, omit_no_match=TRUE)[[1]]
  if(is.null(m)||nrow(m)==0) return(NULL)
  tibble(start=m[,1],end=m[,2])
}

extract_relevant_sections <- function(text){
  
  # FINDINGS to IMPRESSION block extraction strict tolerant colon
  findings_pattern <- "FINDINGS\\s*:\\s*(.*?)(?=IMPRESSION\\s*:)"
  
  findings <- stringr::str_extract(text, findings_pattern)
  
  # IMPRESSION section till end
  impression_pattern <- "IMPRESSION\\s*:\\s*(.*)$"
  
  impression <- stringr::str_extract(text, impression_pattern)
  
  clean <- paste(findings, impression, sep="\n")
  
  # if either was not found fallback to full text (rare safety)
  if(is.na(clean)) return(text)
  
  clean
}



build_phrase_regex <- function(phrases){
  esc <- str_replace_all(phrases,"([\\^$.|?*+()\\[\\]{}\\\\])","\\\\\\1")
  esc <- str_replace_all(esc,"\\s+","\\\\s+")
  paste0("(?i)\\b(",paste(esc,collapse="|"),")(\\b|[,.;:])")
}

split_sentences_with_offsets <- function(text){
  bounds <- stringi::stri_locate_all_boundaries(text,type="sentence")[[1]]
  if(nrow(bounds)==0||any(is.na(bounds))){
    return(tibble(sent_id=1L,sent_text=text,start_char=1L,end_char=nchar(text)))
  }
  tibble(
    sent_id=seq_len(nrow(bounds)),
    sent_text=str_sub(text,bounds[,1],bounds[,2]),
    start_char=bounds[,1],
    end_char=bounds[,2]
  )
}

resolve_overlaps <- function(df_spans){
  df_spans %>%
    arrange(doc_id,sent_id,desc(end_char-start_char+1)) %>%
    group_by(doc_id,sent_id) %>%
    reframe({
      kept<-logical(n());taken<-matrix(FALSE,nrow=0,ncol=2)
      for(i in seq_len(n())){
        s<-start_char[i]; e<-end_char[i]
        overlap<-if(nrow(taken)>0) any(!(e<taken[,1]|s>taken[,2])) else FALSE
        if(!overlap){kept[i]<-TRUE;taken<-rbind(taken,c(s,e))}
      }
      cur_data()[kept,,drop=FALSE]
    })%>%ungroup()
}

scopes_from_triggers <- function(text,trig_df,term_df){
  if(is.null(trig_df)) return(tibble(trig_start=integer(),trig_end=integer(),scope_start=integer(),scope_end=integer()))
  if(is.null(term_df)) term_df<-tibble(start=Inf,end=Inf)
  out<-vector("list",nrow(trig_df))
  for(i in seq_len(nrow(trig_df))){
    t_end<-trig_df$end[i]
    after<-which(term_df$start>t_end)
    scope_end<-if(length(after)>0) term_df$start[min(after)]-1L else nchar(text)
    out[[i]]<-tibble(trig_start=trig_df$start[i],trig_end=t_end,scope_start=t_end+1L,scope_end=scope_end)
  }
  bind_rows(out)
}

assign_context_certainty_sentence <- function(sent_text, ents_sentence_df){
  
  if(nrow(ents_sentence_df)==0) return(ents_sentence_df)
  obs_idx <- which(ents_sentence_df$entity_type=="Observation")
  if(length(obs_idx)==0) return(ents_sentence_df)
  
  pre_neg  <- locate_all(sent_text, rx_ctx$pre_neg)
  post_neg <- locate_all(sent_text, rx_ctx$post_neg)
  pre_unc  <- locate_all(sent_text, rx_ctx$pre_unc)
  post_unc <- locate_all(sent_text, rx_ctx$post_unc)
  term     <- locate_all(sent_text, rx_ctx$term)
  pseudo   <- locate_all(sent_text, rx_ctx$pseudo)
  
  # remove pseudo-negations from neg triggers
  if(!is.null(pre_neg) && !is.null(pseudo)){
    keep <- !sapply(seq_len(nrow(pre_neg)), function(i) any(!(pre_neg$end[i]<pseudo$start | pre_neg$start[i]>pseudo$end)))
    pre_neg <- pre_neg[keep,,drop=FALSE]
    if(nrow(pre_neg)==0) pre_neg <- NULL
  }
  
  # scopes
  neg_scopes <- scopes_from_triggers(sent_text, pre_neg, term)
  unc_scopes <- scopes_from_triggers(sent_text, pre_unc, term)
  
  ents_sentence_df$certainty <- ifelse(ents_sentence_df$entity_type=="Observation","present",NA_character_)
  
  in_any_scope <- function(start_pos, scopes_df){
    if(is.null(scopes_df) || nrow(scopes_df)==0) return(FALSE)
    any(start_pos >= scopes_df$scope_start & start_pos <= scopes_df$scope_end)
  }
  
  # neg scopes
  if(!is.null(neg_scopes) && nrow(neg_scopes)>0){
    for(k in obs_idx){
      s <- ents_sentence_df$start_char[k]
      if(in_any_scope(s, neg_scopes)) ents_sentence_df$certainty[k] <- "absent"
    }
  }
  
  # uncertainty scopes
  if(!is.null(unc_scopes) && nrow(unc_scopes)>0){
    for(k in obs_idx){
      if(ents_sentence_df$certainty[k] == "absent") next
      s <- ents_sentence_df$start_char[k]
      if(in_any_scope(s, unc_scopes)) ents_sentence_df$certainty[k] <- "uncertain"
    }
  }
  
  # post-neg triggers
  if(!is.null(post_neg) && nrow(post_neg)>0){
    for(i in seq_len(nrow(post_neg))){
      t_start <- post_neg$start[i]
      punct <- stringi::stri_locate_all_regex(sent_text,"[,;.:]")[[1]]
      left_bound <- 1L
      if(!is.null(punct) && nrow(punct)>0){
        lefts <- punct[,1][punct[,1]<t_start]
        if(length(lefts)>0) left_bound <- max(lefts)+1L
      }
      for(k in obs_idx){
        if(ents_sentence_df$end_char[k] < t_start && ents_sentence_df$start_char[k] >= left_bound)
          ents_sentence_df$certainty[k] <- "absent"
      }
    }
  }
  
  # post uncertainty triggers
  if(!is.null(post_unc) && nrow(post_unc)>0){
    for(i in seq_len(nrow(post_unc))){
      t_start <- post_unc$start[i]
      punct <- stringi::stri_locate_all_regex(sent_text,"[,;.:]")[[1]]
      left_bound <- 1L
      if(!is.null(punct) && nrow(punct)>0){
        lefts <- punct[,1][punct[,1]<t_start]
        if(length(lefts)>0) left_bound <- max(lefts)+1L
      }
      for(k in obs_idx){
        if(ents_sentence_df$certainty[k]!="absent" &&
           ents_sentence_df$end_char[k]<t_start &&
           ents_sentence_df$start_char[k]>=left_bound)
          ents_sentence_df$certainty[k] <- "uncertain"
      }
    }
  }
  
  # SAFETY FALLBACK
  if(!"certainty" %in% colnames(ents_sentence_df)){
    ents_sentence_df$certainty <- ifelse(ents_sentence_df$entity_type=="Observation","present",NA_character_)
  }
  
  return(ents_sentence_df)
}

extract_entities <- function(texts, doc_ids=NULL,
                             anatomy_dict=anatomy_terms,
                             observation_dict=observation_terms){
  
  if(is.null(doc_ids)) doc_ids <- seq_along(texts)
  
  rx_anat <- build_phrase_regex(unique(anatomy_dict))
  rx_obs  <- build_phrase_regex(unique(observation_dict))
  
  results <- vector("list", length(texts))
  
  for(i in seq_along(texts)){
    raw <- texts[i]
    doc_id <- doc_ids[i]
    focus <- extract_relevant_sections(raw)
    sents <- split_sentences_with_offsets(focus)
    
    spans <- purrr::pmap_dfr(sents, function(sent_id, sent_text, start_char, end_char){
      anat_hits <- stringi::stri_locate_all_regex(sent_text, rx_anat, omit_no_match=TRUE)[[1]]
      obs_hits  <- stringi::stri_locate_all_regex(sent_text, rx_obs,  omit_no_match=TRUE)[[1]]
      
      anat_df <- if(!is.null(anat_hits) && nrow(anat_hits)>0){
        tibble(entity_text=str_sub(sent_text, anat_hits[,1], anat_hits[,2]),
               start_char=start_char+anat_hits[,1]-1L,
               end_char=start_char+anat_hits[,2]-1L,
               sent_id=sent_id, entity_type="Anatomy", sent_text_full=sent_text)
      } else tibble()
      
      obs_df <- if(!is.null(obs_hits) && nrow(obs_hits)>0){
        tibble(entity_text=str_sub(sent_text, obs_hits[,1], obs_hits[,2]),
               start_char=start_char+obs_hits[,1]-1L,
               end_char=start_char+obs_hits[,2]-1L,
               sent_id=sent_id, entity_type="Observation", sent_text_full=sent_text)
      } else tibble()
      
      bind_rows(anat_df, obs_df) %>% mutate(doc_id=doc_id)
    })
    
    if(nrow(spans)==0){ results[[i]] <- tibble(); next }
    
    spans <- spans %>%
      relocate(doc_id, sent_id) %>%
      resolve_overlaps() %>%
      group_by(doc_id, sent_id) %>%
      group_modify(~ assign_context_certainty_sentence(.x$sent_text_full[1], .x)) %>%
      ungroup()
    
    # enforce certainty column presence ALWAYS
    if(!"certainty" %in% colnames(spans)){
      spans$certainty <- NA_character_
    }
    
    spans <- spans %>%
      select(doc_id, sent_id, entity_text, entity_type, certainty, start_char, end_char, sent_text_full) %>%
      arrange(start_char)
    
    
    results[[i]] <- spans
  }
  results <- lapply(results, function(x){
    if(nrow(x)==0){
      tibble(doc_id=character(), sent_id=integer(), entity_text=character(),
             entity_type=character(), certainty=character(), start_char=integer(),
             end_char=integer(), sent_text_full=character())
    } else x
  })
  
  ents <- bind_rows(results)
  
  if(nrow(ents)==0) return(ents)
  
  ents <- ents %>%
    select(-sent_text_full) %>%
    arrange(doc_id, sent_id, start_char) %>%
    mutate(entity_id=row_number(), .before=1L)
  
  ents
}
