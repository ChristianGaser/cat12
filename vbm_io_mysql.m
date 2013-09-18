function varargout=vbm_io_mysql(opt,scans)
% ______________________________________________________________________
% WARNING: THIS IS A PRIVATE FUNCTION THAT IS STILL IN DEVELOPMENT! 
%          IT ONLY WORKS ON OUR SERVERS, AND WILL MAYBE BECOME PART OF
%          THE VBM TOOLBOX IN THE FUTURE.
% ______________________________________________________________________
%
% function varargin=vbm_io_mysql(opt[,scans])
%
% % mySQL datenbank im/export:
%   linux> mysqldump database -u USER -p [tables] > database.sql
%   linux> mysql -p dbname -u USER -p < dumpfile.sql
%
% % Start des Scripts um 02.55 jede Nacht
% #!/bin/sh
% # Backup Script
% # Author: Jerome Griessmeier
% # Version: 0.2
% #
% # This Shell Script backup your database
% # For automating use a cronjob
% 
% #
% # Pfade setzen/ Setting path variables
% #
% MYSQL_DUMP=/usr/bin/mysqldump
% BACKUP_DIR=/pfad/zu/backup_verzeichnis
% TAR=/bin/tar
% RM=/bin/rm
% DB_NAME=DB_NAME
% DB_USER=DB_USER
% DB_PASS=DB_PASS
% AKT_DATUM=`date +%Y%m%d%H%M`
% 
% #
% # mysql dump erzeugen / create mysql dump
% #
% $MYSQL_DUMP $DB_NAME -u $DB_USER --password=$DB_PASS >
% $BACKUP_DIR/$AKT_DATUM.backup.sql
% 
% #
% # mysql dump komprimieren / Compress data
% #
% cd $BACKUP_DIR
% $TAR -cvzf $AKT_DATUM.backup.sql.tgz $AKT_DATUM.backup.sql
% 
% #
% # aufraeumen / clean up
% #
% $RM $AKT_DATUM.backup.sql 
%   linux> crontab -e 55 2 * * * root /backup/backup.sh >> /dev/null 2>&1 
  
  % check if database toolbox exist
  if ~exist('opt','var'), opt=struct(); end
  if exist('scans','var')
    action='insert';
  else
    action='import';
  end

  def.path   = '/usr/local/mysql-5.5.28-osx10.6-x86_64/';
  def.user   = 'root';
  def.pwd    = 'sql4NISALS';
  def.db     = 'vbmDB';
  def.tname  = 'Scan';
  def.tid    = 'ScanID';
  def.url    = '';
  def.DBTB   = 0; 
  opt        = checkinopt(opt,def);
  
  if isempty(opt.url)
    opt.url2   = sprintf('jdbc:mysql://%s/information_schema','localhost:3306');
  else
    opt.url2   = sprintf('jdbc:mysql://%s/information_schema',opt.url);
  end
  
  if opt.DBTB
    try
      conn = database('information_schema',opt.user,opt.pwd,'com.mysql.jdbc.Driver',opt.url2);
      if isempty(conn.url), error('NoDriver'); end
    catch  %#ok<CTCH>
      opt.DBTB = 0;
    end
  end
  if ~opt.DBTB
    if ~isempty(opt.url), url = ['-x ' opt.url]; else url = ''; end
    mysqlcall = sprintf('PATH=$PATH:%s%sbin; mysql %s -u %s -p%s -e',opt.path,filesep,url,opt.user,opt.pwd);
  end
    
  switch action
    case 'insert'
      % create database
      if opt.DBTB
        exec(conn,sprintf('CREATE DATABASE IF NOT EXISTS %s;',opt.db)); %USE %s;,opt.db
        conn = database(opt.db,opt.user,opt.pwd,'com.mysql.jdbc.Driver',opt.url2);
      else
        [SS,SR]=system(sprintf('%s "%s";',mysqlcall,sprintf('CREATE DATABASE IF NOT EXISTS %s;',opt.db)));
      end

      % delete previous table
      if scans(1).(opt.tid)==1
        if opt.DBTB
          e=exec(conn,sprintf('DROP TABLE %s;',opt.tname)); %#ok<NASGU>
        else
          [SS,SR]=system(sprintf('%s "%s";',mysqlcall,sprintf('USE %s; DROP TABLE %s;',opt.db,opt.tname)));
        end
      end
      
      % create new table
      % fields of the table
      fields  = fieldnames(scans);
      if ischar(scans(1).(opt.tid)),  tfields = sprintf('%s CHAR(%d) PRIMARY KEY',opt.tid,255);
      else                            tfields = sprintf('%s INT(%d)  PRIMARY KEY',opt.tid,255);
      end
      if numel(fields)>0
        sfields = setdiff(fields,opt.tid);
        for f=1:numel(sfields)
          [tmp,mysqltype] = matlab2sqltype(scans(1).(sfields{f}));
          tfields = sprintf('%s,%s %s',tfields,sfields{f},mysqltype);
        end
      else
        error('MATLAB:vbm_io_mysql_export:scans contain no fields for export!\n')
      end
      % group is not a valid field!!!
      if opt.DBTB 
        e=exec(conn,sprintf('CREATE TABLE %s (%s);',opt.tname,tfields)); %#ok<NASGU> IF NOT EXISTS 
      else
        [SS,SR]=system(sprintf('%s "%s";',mysqlcall,sprintf('USE %s; CREATE TABLE %s (%s);',opt.db,opt.tname,tfields)));
      end
      
      
      % insert/update data
      fields  = fieldnames(scans);
      for si=1:numel(scans)
        trynr=0; maxtrynr=10; % only for testing
        while trynr<maxtrynr
          if opt.DBTB
            try
              fastinsert(conn,opt.tname,fields,shiftdim(struct2cell(scans(si)),2)');
            catch err 
              if strcmp(err.identifier,'database:database:insertError') && trynr<maxtrynr
                field = textscan(err.message,'%s','delimiter',''''); field=field{1}{2}; 
                [tmp,tfield] = matlab2sqltype(scans(1).(field));
                e=exec(conn,sprintf('ALTER TABLE %s ADD COLUMN %s %s;',opt.tname,field,tfield)); %#ok<NASGU>
                trynr=trynr+1;
              elseif ~isempty(strfind(err.message,'Duplicate entry')) && trynr<maxtrynr
                update(conn,opt.tname,fields,shiftdim(struct2cell(scans),2)',...
                  sprintf('where %s = %d;',opt.tid,scans(si).(opt.tid)));  
                trynr=inf;
              else
                trynr=inf;
                %error('MATLAB:vbm_io_mysql:fastinsert','%s',err.message)
              end
            end
          else
            tflds=''; tvals='';
            for f=1:numel(fields)
              tflds=sprintf('%s,%s',tflds,fields{f});
              if isempty(scans(si).(fields{f}))
                tvals     = sprintf('%s,%s',tvals,'NULL');
              else
                switch class(scans(si).(fields{f}))
                  case {'logical','int','int8','int16','int32','int64','uint8','uint16','uint32','uint64'}
                    tvals = sprintf('%s,%0d',tvals,scans(si).(fields{f}));
                  case {'single','double'}
                    if isnan(scans(si).(fields{f}))
                      tvals = sprintf('%s,%0.6f',tvals,0);
                    else
                      tvals = sprintf('%s,%0.6f',tvals,scans(si).(fields{f}));
                    end
                  case 'char'
                    tvals   = sprintf('%s,%s',tvals,sprintf('''%s''',scans(si).(fields{f})));
                end  
              end
            end
            tflds(1)=[];tvals(1)=[];

            [SS,SR]=system(sprintf('%s "%s";',mysqlcall,sprintf('USE %s; REPLACE INTO %s (%s) VALUES (%s);',opt.db,opt.tname,tflds,tvals)));
             
            if strfind(SR,'ERROR'); 
              if strcmp(SR(1:10),'ERROR 1054')
                field = textscan(SR,'%s','delimiter',''''); field=field{1}{2};
                [tmp,afield] = matlab2sqltype(scans(1).(field));
                [SS,SR]=system(sprintf('%s "%s";',mysqlcall,sprintf('USE %s; ALTER TABLE %s ADD COLUMN %s %s;',opt.db,opt.tname,field,afield))); 
                if strfind(SR,'ERROR'); 
                  error('MATLAB:vbm_io_mysql',SR);
                end
                
              else
                error('MATLAB:vbm_io_mysql',SR);
                 trynr=inf;
              end
            end
            trynr=trynr+1;
          end
        end
      end
    case 'update'
      if opt.DBTB
       update(conn,def.tname,fieldnames(scans),shiftdim(struct2cell(scans),2)',...
          sprintf('where %s = %d;',opt.tid,cell2double(shiftdim(struct2cell(scans),2))));
      else
        % erstmal unwichtig
      end
    case 'import'
      % add fields and scans
      if opt.DBTB
        oldfields = eval(['{' columnnames(fetch(exec(conn,sprintf('SELECT * FROM %s;',opt.tname)))) '}']);
        varargout{1} = cell2struct(fetch(conn,sprintf('SELECT * FROM %s;',opt.tname))',oldfields);
      else
        % unwichtig
      end
    otherwise
      error('MATLAB:vbm_io_mysql:unkown action');
  end
end

function [var,sqltype] = matlab2sqltype(var)
  switch class(var)
    case {'char','cell'}
      var = char(var)'; var = var(:)';
    otherwise
      if isempty(var)
        var = 'NULL';
      else
        var = var(1);
      end
  end
  
  switch class(var)
    case 'logical',       sqltype = 'BOOL';
    case 'uint8',         sqltype = 'TINYINT UNSIGNED';
    case 'uint16',        sqltype = 'SMALLINT UNSIGNED';
    case 'uint32',        sqltype = 'INT UNSIGNED';
    case 'uint64',        sqltype = 'BIGINT UNSIGNED';
    case 'int8',          sqltype = 'TINYINT';
    case 'int16',         sqltype = 'SMALLINT';
    case 'int32',         sqltype = 'INT';
    case 'int64',         sqltype = 'BIGINT';
    case 'single',        sqltype = 'FLOAT';
    case 'double',        sqltype = 'DOUBLE';
    case 'char',          sqltype = 'TEXT';
    case 'cell',          sqltype = 'TEXT';
    otherwise,            error('MATLAB:vbm_io_mysql:matlab2sqltype','ERROR: unsupported type ''%s''',class(var)); 
  end

end

%{
  opt.initpath = sprintf('export set MYSQL_HOME=%s; export set PATH=$PATH:$MYSQL_HOME/bin',opt.path);
  opt.login    = sprintf('mysql -u%s -p%s',opt.user,opt.pwd);
  opt.logout   = 'exit; exit;';
  
  
function sqlc = vbm_io_mysql_export(opt,scans)
  
  % fields of the table
  fields  = fieldnames(scans);
  tfields = 'scanID NOT NULL';
  if numel(fields)>0
    for f=1:numel(fields)
      [tmp,mysqltype] = fmatlab2sqltype(scans(1).(fields{f}));
      tfields = sprintf('%s,%s %s',tfields,scans(1).(fields{f}),mysqltype);
    end
  else
    error('MATLAB:vbm_io_mysql_export:scans contain no fields for export!\n')
  end

  % create database if not exist, use it, and create a table for the scan if not exist
  sqlc = [ ...
    sprintf('CREATE DATABASE IF NOT EXISTS %s;',opt.db) ...      
    sprintf('USE %s;',opt.db) ...                                                   
    sprintf('CREATE TABLE IF NOT EXISTS %s (%s %s);',opt.scans,tfields) ...   
    ];
  
  % export entries
  for s=1:numel(scans)
    sqlc = [ sqlc ...
      sprintf('INSERT INTO %s (%s) VALUES (%s);',tval,tfields{f},tscans(s).(fields{f}) %#ok<AGROW>
      ];
  end
  
end

function subs = vbm_io_mysql_import()
  % read entry
end
%}