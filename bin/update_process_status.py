#!/usr/bin/env python

import argparse
from sqlalchemy import create_engine, update
from sqlalchemy.ext.automap import automap_base
from sqlalchemy.orm import Session

parser = argparse.ArgumentParser(description="Update process status for database entry mode")
parser.add_argument('-d','--database',dest='database', help="database name")
parser.add_argument('-H','--hostname',dest='hostname')
parser.add_argument('-u','--user',dest='dbuser')
parser.add_argument('-p','--password',dest='mariadbpassword',)
parser.add_argument('--status', dest='status', help='new status to be set')
parser.add_argument('--ids', dest='ids', help='comma separated list of process ids affected')

def main():
    args = parser.parse_args()

    status = args.status
    ids = [int(i.strip()) for i in args.ids.split(",")]
    
    engine = create_engine("mysql+mysqlconnector://%s:%s@%s:3306/%s" %
             (args.dbuser, args.mariadbpassword, args.hostname,args.database))
    Base = automap_base()
    Base.prepare(autoload_with=engine)
    
    Process = Base.classes.process

    with Session(engine) as session:
        stmt = update(Process).where(Process.id.in_(ids)).values(status=status)
        session.execute(stmt)
        session.commit()

if __name__ == "__main__":
    main()

