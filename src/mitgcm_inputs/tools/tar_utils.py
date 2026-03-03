import grp
import os
import pwd
import tarfile
from functools import lru_cache
from time import time


@lru_cache(maxsize=1)
def get_current_user_description():
    """
    Retrieves information about the current user.

    This function gathers details about the current user by accessing
    system-level information such as user and group IDs and resolving their
    respective names via system libraries. It is cached for efficiency using
    an LRU cache.

    Returns:
        A dictionary containing the following keys:
            - uid (int): The user ID of the current user.
            - gid (int): The group ID of the current user.
            - uname (str): The username of the current user.
            - gname (str): The group name of the current user's group.
    """
    uid = os.getuid()
    gid = os.getgid()
    return {
        "uid": uid,
        "gid": gid,
        "uname": pwd.getpwuid(uid).pw_name,
        "gname": grp.getgrgid(gid).gr_name,
    }


def set_tar_file_ownerships(tar_pointer: tarfile.TarInfo):
    """
    Sets the ownership and modification time for a tar file member.

    This function updates the ownership attributes (`uid`, `gid`, `uname`,
    `gname`) and the modification time (`mtime`) of a tarfile.TarInfo object
    using details obtained from the current user's description.

    Args:
        tar_pointer: A tarfile.TarInfo object representing the tar file member
        whose ownership details and modification time will be updated.
    """
    user_description = get_current_user_description()

    tar_pointer.uid = user_description["uid"]
    tar_pointer.gid = user_description["gid"]
    tar_pointer.uname = user_description["uname"]
    tar_pointer.gname = user_description["gname"]
    tar_pointer.mtime = int(time())
