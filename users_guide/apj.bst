%$$$ apj.bst $$$
% This file is an output file created from astrobib.bst.  Please edit
% that file if you have any changes.
%
% This style file is for BibTeX version 0.99c, LaTeX version 2.09
% Place it in a file called apj.bst in the BibTeX search path.
% (Placing it in the same directory as the LaTeX document should also work.)
% You need to include astrobib.sty as a \documentstyle option.









% 
% ApJ BibTeX style -- ApJ style with annotations.
% Articles are cited in text as in chicago.bst and chicagoa.bst, described
% below.
%
% In the references section, articles are cited as described in

% Abt 1990 (ApJ 357, 1), with some of the ambiguities cleared up
% by looking at examples from recent ApJ's (not AJ's, which were formatted
% wrong by the publishing house).

%
% Please report bugs to ferguson@stsci.edu
% =====================================================================
% Revision history:
% v1  17 May 1993       H. Ferguson     Adapted from chicagoa.bst
% v2  25 May 1993       H. Ferguson     Fixed long name bug
% v3  06 Oct 1994       E. Bradford     Many corrections to ApJ style
%     10 Oct 1994       E. Bradford     Merged versions to one master file.
%     20 Oct 1994       E. Bradford     Many updates for non-ApJ journals.
%
% ====================================================
% Original documentation from chicagoa.bst
%
% This "chicago" style is based on newapa.bst (American Psych. Assoc.)
% found at ymir.claremont.edu. Annotation support added 07/09/92.
% Modifications are loosely based on Chicago Manual of Style.
%
%   Citation format: (author-last-name year)
%             (author-last-name and author-last-name year)
%             (author-last-name, author-last-name, and author-last-name year)
%             (author-last-name et al. year)
%             (author-last-name)
%             author-last-name (year)
%             (author-last-name and author-last-name)
%             (author-last-name et al.)
%             (year) or (year,year)
%             year or year,year
%
%   Reference list ordering: alphabetical by author or whatever passes
%    for author in the absence of one.
%
% This BibTeX style has support for abbreviated author lists and for
%    year-only citations.  This is done by having the citations
%    actually look like
%
%    \citeauthoryear{full-author-info}{abbrev-author-info}{year}
%
% The LaTeX style has to have the following (or similar)
%
%     \let\@internalcite\cite
%     \def\fullcite{\def\citeauthoryear##1##2##3{##1, ##3}\@internalcite}
%     \def\fullciteA{\def\citeauthoryear##1##2##3{##1}\@internalcite}
%     \def\shortcite{\def\citeauthoryear##1##2##3{##2, ##3}\@internalcite}
%     \def\shortciteA{\def\citeauthoryear##1##2##3{##2}\@internalcite}
%     \def\citeyear{\def\citeauthoryear##1##2##3{##3}\@internalcite}
%
% These TeX macro definitions are found in chicago.sty. Additional
% commands to manipulate different components of a citation can be defined
% so that, for example, you can list author's names without parentheses
% if using a citation as a noun or object in a sentence.
%
% This file was originally copied from newapa.bst at ymir.claremont.edu.
%
% Features of aasa.bst:
% =========================
%
% - supports an "annotation" field for annotated bibliographies.
% - full names used in citations, but abbreviated citations are available
%   (see above)
% - if an entry has a "month", then the month and year are also printed
%   as part of that reference item.
% - book, inbook and manual use "location: publisher" (or organization)
%   for address and publisher. All other types list publishers separately.
% - "pp." are used to identify page numbers for all entry types except
%   articles.
% - organization is used as a citation label if neither author nor editor
%   is present (for manuals).
% - "et al." is used for long author and editor lists.
%
% Modified by Glenn Paulley, University of Waterloo, July, 1992 (mostly
% from Chicago Manual of Style):
%
% Modifications: (from newapa.bst)
% =============
%
%   - added month, year to bib entries if month is present
%   - fixed bug with In proceedings, added necessary comma after title
%   - all conjunctions changed to "and" from "\&"
%   - fixed bug with author labels in my.full.label: "et al." now is
%        generated when "others" is an author name
%   - major modification from Chicago Manual of Style (13th ed.) is that
%     only the first author in a reference appears last name first-
%     additional authors appear as J. Q. Public.
%   - pages are listed as "pp. xx-xx" in all entry types, except
%     article entries. Unnecessary (IMHO) "()" around page numbers
%     were removed, and page numbers now don't end with a period.
%   - created chicago.sty for use with this bibstyle (required).
%   - fixed bugs in FUNCTION {format.vol.num.pages} for missing volume,
%     number, and /or pages. Renamed to format.jour.vol.
%   - fixed bug in formatting booktitles: additional period an error if
%     book has a volume.
%   - fixed bug: editors usually given redundant period before next clause
%     (format.editors.dot) removed.
%   - added label support for organizations, if both author and editor
%     are missing (from alpha.bst). If organization is too long, then
%     the key field is used for abbreviated citations.
%   - In proceedings or books of several volumes, no comma was written
%     between the "Volume x" and the page numbers (this was intentional
%     in newapa.bst). Fixed.
%   - Some journals may not have volumes/numbers, only month/year (eg.
%     IEEE Computer). Fixed bug in article style that assumed volume/number
%     was always present.
%
% Original documentation for newapa.sty:
% =====================================
%
% This version was made by modifying the master file made by
% Oren Patashnik (opbibtex@CS.STANFORD.EDU), and the 'named' BibTeX
% style of Peter F. Patel-Schneider.
%
% Copyright (C) 1985, all rights reserved.
% Copying of this file is authorized only if either
% (1) you make absolutely no changes to your copy, including name, or
% (2) if you do make changes, you name it something other than 'newapa.bst'.
% There are undoubtably bugs in this style.  If you make bug fixes,
% improvements, etc.  please let me know.  My e-mail address is:
%    spencer@cgrg.ohio.state.edu or 71160.3141@compuserve.com
%
% This style was made from 'plain.bst', 'named.bst', and 'apalike.bst',
% with lots of tweaking to make it look like APA style, along with tips
% from Young Ryu and Brian Reiser's modifications of 'apalike.bst'.

ENTRY
  { address
    annotation
    author
    booktitle
    chapter
    edition
    editor
    howpublished
    institution
    journal
    key
    month
    note
    number
    organization
    pages
    publisher
    school
    series
    title
    type
    volume
    year
  }
  {}
  { label l1 l2 l3 extra.label sort.label }

FUNCTION {not}
{   { #0 }
    { #1 }
  if$
}

FUNCTION {and}
{   'skip$
    { pop$ #0 }
  if$
}

FUNCTION {or}
{   { pop$ #1 }
    'skip$
  if$
}

INTEGERS { output.state before.all mid.sentence after.block 
	   mid.sen.nocomma }

FUNCTION {init.state.consts}
{ #0 'before.all :=
  #1 'mid.sentence :=
  #2 'after.block :=
  #3 'mid.sen.nocomma :=
}

STRINGS { s t u }

FUNCTION {output.nonnull}
{ duplicate$ 's :=
  % If the string to be output starts with punctuation, or if we are
  % at the beginning, don't add any punctuation or space.
  #1 #1 substring$ % First character in the new string just added.
  duplicate$ % This will be used to check for an open parenthesis.
  duplicate$ "," = swap$
  duplicate$ "." = swap$
  duplicate$ ";" = swap$
  duplicate$ ":" = swap$
  duplicate$ "!" = swap$
  "?" = or or or or or
  output.state before.all = or
    'pop$ % Eliminate first character which would have been used for paren chk.
    {
     "(" =
      output.state mid.sen.nocomma = or
	{ " " * }
	{ output.state after.block =
	  {  add.period$ write$
	     newline$
	     "\newblock "
	  }
	  {  output.state mid.sentence =
	     { ", " * }
	     { "Internal error: unrecognized state in " cite$ * warning$ }
	     if$
	  }
	  if$
	}
	if$
    }
    if$
  write$
  mid.sentence 'output.state :=
  s
}

% Use a colon to separate output. Used only for address/publisher
% combination in book/inbook types.
%
%FUNCTION {output.nonnull.colon}
%{ 's :=
%  output.state mid.sentence =
%  { ": " * write$ }
%  {  output.state after.block =
%     {  add.period$ write$
%        newline$
%        "\newblock " write$
%     }
%     {  output.state before.all =
%        'write$
%        {  output.state mid.sen.nocomma =
%           {  " " * write$ }
%           {  add.period$ " " * write$ }
%           if$
%        }
%        if$
%     }
%     if$
%     mid.sentence 'output.state :=
%  }
%  if$
%  s
%}

FUNCTION {output}
{ duplicate$ empty$
    'pop$
    'output.nonnull
  if$
}

FUNCTION {output.check}
{ 't :=
  duplicate$ empty$
    { pop$ "empty " t * " in " * cite$ * warning$ }
    'output.nonnull
  if$
}

%FUNCTION {output.colon}
%{ duplicate$ empty$
%    'pop$
%    'output.nonnull.colon
%  if$
%}

%FUNCTION {output.check.colon}
%{ 't :=
%  duplicate$ empty$
%    { pop$ "empty " t * " in " * cite$ * warning$ }
%    'output.nonnull.colon
%  if$
%}

FUNCTION {output.year}
{ year empty$
     'skip$
     { write$

         " "

         year * extra.label *
       mid.sentence 'output.state :=
     }
  if$
}

FUNCTION {output.year.check}
{ year empty$
     { "empty year in " cite$ * warning$ }
     { output.year }
  if$
}

FUNCTION {output.bibitem}
{ newline$
  "\bibitem[" write$
  l2 l1 "\protect\citeauthoryear{" swap$ * "}{" * swap$ * "}{" * write$
  l3 write$
  "]{" write$
  cite$ write$
  "}" write$
  newline$
  ""
  before.all 'output.state :=
}

FUNCTION {new.block}
{ output.state before.all =
    'skip$
    { after.block 'output.state := }
  if$
}

%FUNCTION {new.sentence}
%{ output.state after.block =
%    'skip$
%    { output.state before.all =
%    'skip$
%    { after.sentence 'output.state := }
%      if$
%    }
%  if$
%}

%FUNCTION {new.block.checka}
%{ empty$
%    'skip$
%    'new.block
%  if$
%}

%FUNCTION {new.block.checkb}
%{ empty$
%  swap$ empty$
%  and
%    'skip$
%    'new.block
%  if$
%}

%FUNCTION {new.sentence.checka}
%{ empty$
%    'skip$
%    'new.sentence
%  if$
%}

%FUNCTION {new.sentence.checkb}
%{ empty$
%  swap$ empty$
%  and
%    'skip$
%    'new.sentence
%  if$
%}

FUNCTION {field.or.null}
{ duplicate$ empty$
    { pop$ "" }
    'skip$
  if$
}

%FUNCTION {underline}
%{ duplicate$ empty$
%  { pop$ "" }
%  { "\underline{" swap$ * "}" * }
%  if$
%}

%
% Emphasize the top string on the stack.
%
% FUNCTION {emphasize}
% { duplicate$ empty$
%    { pop$ "" }
%    { "{\em " swap$ * "}" * }
%  if$
%}

%
% Emphasize the top string on the stack, but add a trailing space.
%
% FUNCTION {emphasize.space}
% { duplicate$ empty$
%   { pop$ "" }
%   { "{\em " swap$ * "\/}" * }
% if$
% }

FUNCTION {format.annotation}
{ annotation empty$
	{ "" }
	{ " \begin{quotation}\noindent "
	  annotation
	  * " \end{quotation} " *
	}
      if$
}

FUNCTION {fin.entry}
{
  note output
  write$
  newline$
  format.annotation write$
}

INTEGERS { nameptr namesleft numnames }

%
% Format bibliographical entries with the last name first.
%
FUNCTION {format.names}
{ 's :=   
  #1 'nameptr :=
  s num.names$ 'numnames :=
% Check for "et al." condition.

  numnames #8 >
  { s nameptr 

    "{vv~}{ll}{, f.}{, jj}"

    format.name$

    "," *


    " et~al." * }
  {
    numnames 'namesleft :=
    { namesleft #0 > }
    { s nameptr

      "{vv~}{ll}{, f.}{, jj}"

      format.name$ 't :=
      nameptr #1 >
        { namesleft #1 >
            { ", " * t * }
            { numnames #2 >
	        {", " *}

	        {", " *}

              if$
              t "others" =
                { " et~al." * }
                { 

                  " \& "

                  * t * }
              if$
            }
          if$
        }
        't
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
    while$
  }
  if$
}   
 
% Just last names.
% If more than three authors, use et al.
FUNCTION {my.full.label}
{ 's :=
  #1 'nameptr :=               % nameptr = 1;
  s num.names$ 'numnames :=    % numnames = num.name$(s);
  numnames 'namesleft :=
  namesleft #3 >
  { s nameptr "{vv~}{ll}" format.name$ " et~al." * }
  {
    { namesleft #0 > }
    { s nameptr "{vv~}{ll}" format.name$ 't :=  % get the next name
      nameptr #1 >
        { namesleft #1 >
              { ", " * t * }
              { numnames #2 >
                  { "," * }
                  'skip$
                if$
                t "others" =
                    { " et~al." * }

                    { " \& " * t * }

                  if$
               }
               if$
             }
            't
        if$
        nameptr #1 + 'nameptr :=          % nameptr += 1;
        namesleft #1 - 'namesleft :=      % namesleft =- 1;
    }
    while$
  }
  if$
}

FUNCTION {format.authors}
{ author empty$
    { "" }
    { author format.names }
  if$
}

FUNCTION {format.key}
{ empty$
    { key field.or.null }
    { "" }
  if$
}

FUNCTION {add.editors.indication}
{ editor num.names$ #1 >

    { ", ed" * }
    { ", ed" * }

  if$
}


FUNCTION {format.names.fml}
%
% Format names in "familiar" format, with first initial followed by
% last name. Like format.names, ALL names are formatted.
%
{ 's :=
  #1 'nameptr :=               % nameptr = 1;
  s num.names$ 'numnames :=    % numnames = num.name$(s);
  numnames 'namesleft :=
    { namesleft #0 > }

    { s nameptr "{f.~}{vv~}{ll}{, jj}" format.name$ 't :=

      nameptr #1 >
        { namesleft #1 >
              { ", " * t * }
               { numnames #2 >
                    { "," * }
                    'skip$
                  if$
                  t "others" =
                    { " et~al." * }
                    { 

                      " \& "

                       * t * }
                  if$
                }
               if$
             }
            't
        if$
        nameptr #1 + 'nameptr :=          % nameptr += 1;
        namesleft #1 - 'namesleft :=      % namesleft =- 1;
    }
  while$
}

%
% Format editor names for use in the "in" types: inbook, incollection,
% inproceedings: first initial, then last names. When editors are the
% LABEL for an entry, then format.editor is used which lists editors
% by last name first.
%
FUNCTION {format.editors.fml}
{ editor empty$
    { "" }
    { editor format.names.fml }
  if$
}


%
% Format editor names for use in labels, last names first.
%
FUNCTION {format.editors}
{ editor empty$
    { "" }
    { editor format.names add.editors.indication }
  if$
}

FUNCTION {format.title}
{ title empty$
    { "" }

    { title "t" change.case$ }

  if$
}

FUNCTION {output.publisher}
{

  mid.sen.nocomma 'output.state :=
  publisher empty$
    { address empty$
        { "No publisher or address in " cite$ * warning$ ""}
        { "There's an address but no publisher in " cite$ * warning$
          "(" address * ")" *
        }
      if$
    }
    { address empty$
        { "(" publisher * ")" * }

        { "(" address * ": " * publisher * ")" * }

      if$
    }
  if$
  output

}

FUNCTION {tie.or.space.connect}
{ duplicate$ text.length$ #3 <
    { "~" }
    { " " }
  if$
  swap$ * *
}

FUNCTION {format.btitle} % Title should be on stack.
{ duplicate$ empty$ edition empty$ or
    'skip$
    { " (" * edition * " ed.)" * }
  if$
}

FUNCTION {format.ed.booktitle} % The title should be on the stack.
{ duplicate$ empty$
  { "no book title in " cite$ * warning$ "" pop$ }
  { editor empty$
    author empty$ or  % Empty author means editor already given.
      'format.btitle

      { format.btitle ", ed. " * format.editors.fml * }

      if$
  }
  if$
}

FUNCTION {add.no}
{ " No." number tie.or.space.connect *
}

FUNCTION {add.vol}
{ ", Vol." volume tie.or.space.connect *
}

FUNCTION {format.full.book.spec} % The title should be on the stack.
{ series empty$
    { format.ed.booktitle
      volume empty$
	{ number empty$
	    'skip$
	    { " there's a number but no series in " cite$ * warning$
	      add.no }
	    if$
	}
	{ add.vol
	  number empty$
	    'skip$
	    {"Both volume and number fields in " * cite$ * warning$ }
	    if$
	}
	if$
    }
    { volume empty$
	{ format.ed.booktitle ", " * series *
	  number empty$
	    'skip$
	    'add.no
	    if$
	}
	{ series add.vol ", " * swap$ format.ed.booktitle *
	  number empty$
	    'skip$
	    {"Both volume and number fields in " * cite$ * warning$ }
	    if$
	}
        if$
    }
    if$
}

INTEGERS { multiresult }

%FUNCTION {multi.page.check}
%{ 't :=
%  #0 'multiresult :=
%    { multiresult not
%      t empty$ not
%      and
%    }
%    { t #1 #1 substring$
%      duplicate$ "-" =
%      swap$ duplicate$ "," =
%      swap$ "+" =
%      or or
%    { #1 'multiresult := }
%    { t #2 global.max$ substring$ 't := }
%      if$
%    }
%  while$
%  multiresult
%}

%Return only the first page from page range (e.g. "1" from "1-20").
FUNCTION {format.pages.short}
{ pages 't :=
  ""
  { t empty$ not
    t #1 #1 substring$ "-" = not
    and
    t #1 #1 substring$ "+" = not
    and
  }
    { t #1 #1 substring$ *
      t #2 global.max$ substring$ 't :=
    }
  while$
}

% Format pages for a book citation.
FUNCTION {format.pages.p}
{
  pages empty$
    { "" }
    { "p." format.pages.short tie.or.space.connect }
  if$
}

FUNCTION {format.pages.book}
{
  pages empty$
    { "" }

    'format.pages.short

  if$
}

% By Young (and Spencer)
% GNP - fixed bugs with missing volume, number, and/or pages
%
% Format journal, volume, number, pages for article types.
%
FUNCTION {output.jour.vol}
{ journal empty$
    { "no journal in " cite$ * warning$
      "" }
    { journal output }
  if$
  volume empty$
    'skip$
    { Volume output }
  if$
  pages empty$
    'skip$
    { 
      duplicate$ empty$ % For no volume, add "p." so page different from vol.
	{ pop$ ", " * format.pages.p * }
	{ ", " * format.pages.short * }
      if$
      
    }
  if$
}

FUNCTION {format.chapter.pages}
{ chapter empty$
    'format.pages.book
    { type empty$
        { "Ch." }
        { type "t" change.case$ }
      if$
      chapter tie.or.space.connect
      pages empty$
        'skip$
        { ", " * format.pages.book * }
      if$
    }
  if$
}

FUNCTION {format.thesis.type}
{ type empty$
    'skip$
    { pop$
      type "t" change.case$
    }
  if$
}

FUNCTION {format.tr.number}
{ type empty$
    { "Technical Report" }
    'type
    if$
  number empty$
    { "t" change.case$ }
    { number tie.or.space.connect }
    if$
}

%FUNCTION {format.crossref}
%{ "in \citeN{" crossref * "}" *
%}

%FUNCTION {format.book.crossref}
%{ volume empty$
%    { "empty volume in " cite$ * "'s crossref of " * crossref * warning$
%      "in "
%    }
%    { "Vol." volume tie.or.space.connect
%      " of " *
%    }
%  if$
%  " \citeN{" * crossref * "}" *
%}

FUNCTION {article}
{ output.bibitem
  format.authors output
  author format.key output          % added

  add.period$

  output.year.check               
  output.jour.vol
  fin.entry
}

FUNCTION {book}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { editor empty$
        { format.authors output.nonnull }
	{ "Can't use both author and editor fields in " * cite$ * warning$ }
        if$
    }
  if$

  add.period$

  output.year.check       % added
  title format.full.book.spec output
  output.publisher
  format.pages.book output
  fin.entry
}

FUNCTION {booklet}
{ output.bibitem
  format.authors output
  author format.key output

  author empty$ key empty$ and
    'skip$
    'add.period$
  if$
  mid.sen.nocomma 'output.state :=

  output.year
  format.title output
  howpublished output
  address output
  fin.entry
}

FUNCTION {inbook}
{ output.bibitem
  author empty$
    { format.editors "author and editor" output.check }
    { format.authors output.nonnull }
    if$

  add.period$

  output.year.check                 % added
  title format.full.book.spec 

  output
  output.publisher
  format.chapter.pages "chapter and pages" output.check
  fin.entry
}

FUNCTION {incollection}
{ output.bibitem
  format.authors "author" output.check
  author format.key output       % added

  add.period$

  output.year.check              % added
  "in " booktitle format.full.book.spec * output
  output.publisher
  format.chapter.pages output
  fin.entry
}

FUNCTION {inproceedings}
{ output.bibitem
  format.authors "author" output.check
  author format.key output            % added

  add.period$

  output.year.check                   % added

  "in " booktitle format.full.book.spec * output
  publisher empty$
    { organization output
      address output
    }
    { output.publisher }
    if$
  format.pages.book output
  fin.entry
}

FUNCTION {conference} { inproceedings }

FUNCTION {manual}
{ output.bibitem
  author empty$
    { editor empty$
      { organization format.key output }  % if all else fails, use key
      { format.editors "author and editor" output.check }
      if$
    }
    { format.authors output.nonnull }
    if$

  add.period$
  mid.sen.nocomma 'output.state :=

  output.year
  title format.btitle "title" output.check
  author empty$ not editor empty$ not or
     {organization output}
     'skip$
     if$
  address output
%  address output
%  ":" output
%  organization output
  fin.entry
}

FUNCTION {misc}
{ output.bibitem
  format.authors output
  author format.key output            % added

  author empty$ key empty$ and
    'skip$
    'add.period$
  if$
  mid.sen.nocomma 'output.state :=

  output.year
  format.title output
  howpublished output
  fin.entry
}

FUNCTION {privatecomm}
{ output.bibitem
  format.authors "author" output.check
  author format.key output                              % special for

  add.period$

  output.year.check                                     % apalike
  "private communication" output
  fin.entry
}

FUNCTION {output.thesis}
{ 'u :=			% Thesis type is passed in.
  output.bibitem
  format.authors "author" output.check
  author format.key output          % added

  add.period$

  output.year.check                 % added
  u format.thesis.type output.nonnull
  school "school" output.check
  address output
  fin.entry
}

FUNCTION {mastersthesis}
{ "Master's thesis"
  output.thesis
}

FUNCTION {phdthesis}
{
  "Ph.D. thesis"
  output.thesis
}

FUNCTION {proceedings}
{ output.bibitem
  editor empty$
    { organization output
      organization format.key output }  % gnp - changed from author format.key
    { format.editors output.nonnull }
  if$
% author format.key output             % gnp - removed (should be either
%                                        editor or organization

  add.period$

  output.year.check                    % added (newapa)
  title format.full.book.spec output
  publisher empty$ not      % No need for warning if no pub.
    { output.publisher }
    { editor empty$         % For empty editor, organization was already given.
        'skip$
        { organization output }
        if$
      address output
    }
    if$
  fin.entry
}

FUNCTION {techreport}
{ output.bibitem
  format.authors
  "author" output.check
  author format.key output             % added

  add.period$

  output.year.check                    % added
  format.title "title" output.check
  format.tr.number output.nonnull
  institution "institution" output.check
  address output
  fin.entry
}

FUNCTION {unpublished}
{ output.bibitem
  format.authors "author" output.check
  author format.key output              % added

  add.period$
  mid.sen.nocomma 'output.state :=

  output.year
  fin.entry
}

FUNCTION {default.type} { misc }

MACRO {jan} {"Jan."}

MACRO {feb} {"Feb."}

MACRO {mar} {"Mar."}

MACRO {apr} {"Apr."}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"Aug."}

MACRO {sep} {"Sept."}

MACRO {oct} {"Oct."}

MACRO {nov} {"Nov."}

MACRO {dec} {"Dec."}

READ

FUNCTION {sortify}
{ purify$
  "l" change.case$
}

INTEGERS { len }

FUNCTION {chop.word}
{ 's :=
  'len :=
  s #1 len substring$ =
    { s len #1 + global.max$ substring$ }
    's
  if$
}

% format.lab.names:
%
% determines "short" names for the abbreviated author information.
% "Long" labels are created in calc.label, using the routine my.full.label
% to format author and editor fields.
%
% There are 4 cases for labels.   (n=3 in the example)
% a) one author             Foo
% b) one to n               Foo, Bar and Baz
% c) use of "and others"    Foo, Bar et al.
% d) more than n            Foo et al.
%
FUNCTION {format.lab.names}
{ 's :=
  s num.names$ 'numnames :=
  numnames #2 >    % change number to number of others allowed before
		   % forcing "et al".
    { s #1 "{vv~}{ll}" format.name$ " et~al." * }
    {
      numnames #1 - 'namesleft :=
      #2 'nameptr :=
      s #1 "{vv~}{ll}" format.name$
	{ namesleft #0 > }
	{ nameptr numnames =
	    { s nameptr "{ff }{vv~}{ll}{ jj}" format.name$ "others" =
		{ " et~al." * }
		{

                  " \& "

                  * s nameptr "{vv~}{ll}" format.name$ * }
	      if$
	    }
	    { ", " * s nameptr "{vv~}{ll}" format.name$ * }
	  if$
	  nameptr #1 + 'nameptr :=
	  namesleft #1 - 'namesleft :=
	}
      while$
    }
  if$
}

FUNCTION {author.key.label}
{ author empty$
    { key empty$
          { "no key, author in " cite$ * warning$
            cite$ #1 #3 substring$ }
         'key
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.label}
{ editor empty$
    { key empty$
          { "no key, ed. in " cite$ * warning$
            cite$ #1 #3 substring$ }
          'key
        if$
     }
     { editor format.lab.names }
  if$
}

FUNCTION {author.key.organization.label}
%
% added - gnp. Provide label formatting by organization if author is null.
%
{ author empty$
    { organization empty$
	{ key empty$
	    { "no key, author or org. in " cite$ * warning$
              cite$ #1 #3 substring$ }
	    'key
	  if$
	}
        { organization text.length$ #30 >
          { "org. too long, using key in " cite$ * warning$
            key }
          { organization }
          if$
        }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {editor.key.organization.label}
%
% added - gnp. Provide label formatting by organization if editor is null.
%
{ editor empty$
    { organization empty$
	{ key empty$
	    { "no key, editor or org. in " cite$ * warning$
              cite$ #1 #3 substring$ }
	    'key
	  if$
	}
        { organization text.length$ #30 >
          { "org. too long, using key in " cite$ * warning$
            key }
          { organization }
          if$
        }
      if$
    }
    { editor format.lab.names }
  if$
}

FUNCTION {author.editor.key.label}
{ author empty$
    { editor empty$
          { key empty$
               { "no key, author, or editor in " cite$ * warning$
                 cite$ #1 #3 substring$ }
             'key
           if$
         }
          { editor format.lab.names }
      if$
    }
    { author format.lab.names }
  if$
}

FUNCTION {calc.label}
%
% Changed - GNP. See also author.organization.sort, editor.organization.sort
% Form label for BibTeX entry. The classification of which fields are used
% for which type of entry (book, inbook, etc.) are taken from alpha.bst.
% The change here from newapa is to also include organization as a
% citation label if author or editor is missing.
%
{ type$ "book" =
  type$ "inbook" =
  or
    'author.editor.key.label
    { type$ "proceedings" =
	'editor.key.organization.label
	{ type$ "manual" =
	    'author.key.organization.label
	    'author.key.label
	  if$
	}
      if$
    }
  if$
  duplicate$

  author empty$  % generate the full label citation information.
    { editor empty$
        { organization empty$
	   { key empty$
	      { year empty$
		 { title empty$
		    { "no field for label in " cite$ * warning$
		      "??" }
		    { title }
		    if$
                 }
		 { year }
		 if$
              }
	      { key }
	      if$
           }
           { organization }
           if$
        }
        { editor my.full.label }
        if$
    }
    { author my.full.label }
  if$

  'l1 :=
  'l2 :=
  year field.or.null purify$ #-1 #4 substring$ *
  sortify 'sort.label :=
  year field.or.null purify$ #-1 #4 substring$ 
  'l3 :=
  l1 l2 * l3 * 'label :=
}

FUNCTION {sort.format.names}
{ 's :=
  #1 'nameptr :=
  ""
  s num.names$ 'numnames :=
  numnames 'namesleft :=
    { namesleft #0 > }
    { nameptr #1 >
          { "   " * }
         'skip$
      if$
      s nameptr "{vv{ } }{ll{ }}{  f{ }}{  jj{ }}" format.name$ 't :=
      nameptr numnames = t "others" = and
          { " et al" * }
          { t sortify * }
      if$
      nameptr #1 + 'nameptr :=
      namesleft #1 - 'namesleft :=
    }
  while$
}

FUNCTION {sort.format.title}
{ 't :=
  "A " #2
    "An " #3
      "The " #4 t chop.word
    chop.word
  chop.word
  sortify
  #1 global.max$ substring$
}

FUNCTION {author.sort}
{ author empty$
    { key empty$
         { "to sort, need author or key in " cite$ * warning$
           "" }
         { key sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.sort}
{ editor empty$
    { key empty$
         { "to sort, need editor or key in " cite$ * warning$
           ""
         }
         { key sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {author.editor.sort}
{ author empty$
    { % "missing author in " cite$ * warning$
      editor empty$
         { key empty$
             { "to sort, need author, editor, or key in " cite$ * warning$
               ""
             }
             { key sortify }
           if$
         }
         { editor sort.format.names }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {author.organization.sort}
%
% added - GNP. Stack author or organization for sorting (from alpha.bst).
% Unlike alpha.bst, we need entire names, not abbreviations
%
{ author empty$
    { organization empty$
	{ key empty$
	    { "to sort, need author, org., or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ organization sortify }
      if$
    }
    { author sort.format.names }
  if$
}

FUNCTION {editor.organization.sort}
%
% added - GNP. Stack editor or organization for sorting (from alpha.bst).
% Unlike alpha.bst, we need entire names, not abbreviations
%
{ editor empty$
    { organization empty$
	{ key empty$
	    { "to sort, need editor, org., or key in " cite$ * warning$
	      ""
	    }
	    { key sortify }
	  if$
	}
	{ organization sortify }
      if$
    }
    { editor sort.format.names }
  if$
}

FUNCTION {presort}
%
% Presort creates the bibentry's label via a call to calc.label, and then
% sorts the entries based on entry type. Chicago.bst adds support for
% including organizations as the sort key; the following is stolen from
% alpha.bst.
%
{ calc.label
  label sortify
  "    "
  *
  type$ "book" =
  type$ "inbook" =
  or
    'author.editor.sort
    { type$ "proceedings" =
	'editor.organization.sort
	{ type$ "manual" =
	    'author.organization.sort
	    'author.sort
	  if$
	}
      if$
    }
  if$
  #1 entry.max$ substring$        % added for newapa
  'sort.label :=                  % added for newapa
  sort.label                      % added for newapa
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT             % by label, sort.label, title --- for final label calculation

STRINGS { last.label next.extra }

INTEGERS { last.extra.num }

FUNCTION {initialize.extra.label.stuff}
{ #0 int.to.chr$ 'last.label :=
  "" 'next.extra :=
  #0 'last.extra.num :=
}

FUNCTION {forward.pass}
{ last.label label =
     { last.extra.num #1 + 'last.extra.num :=
       last.extra.num int.to.chr$ 'extra.label :=
     }
     { "a" chr.to.int$ 'last.extra.num :=
       "" 'extra.label :=
       label 'last.label :=
    }
  if$
}

FUNCTION {reverse.pass}
{ next.extra "b" =
    { "a" 'extra.label := }
     'skip$
  if$
  l3 extra.label * "}" * 'l3 :=
  label extra.label * "}" * 'label :=
  extra.label 'next.extra :=
}

EXECUTE {initialize.extra.label.stuff}

ITERATE {forward.pass}

REVERSE {reverse.pass}

FUNCTION {bib.sort.order}
{ sort.label
  "    "
  *
  year field.or.null sortify
  *
  "    "
  *
  title field.or.null
  sort.format.title
  *
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {bib.sort.order}

SORT             % by sort.label, year, title --- giving final bib. order.

FUNCTION {begin.bib}

{ preamble$ empty$
    'skip$
    { preamble$ write$ newline$ }
  if$
  "\begin{thebibliography}{}" write$ newline$
}


EXECUTE {begin.bib}

EXECUTE {init.state.consts}

ITERATE {call.type$}

FUNCTION {end.bib}
{ newline$
  "\end{thebibliography}" write$ newline$
}

EXECUTE {end.bib}